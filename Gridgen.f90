Module GridGen

    Type Stencil
        Real*8  :: e , w , n , s , p
        Real*8  :: se, sw, ne, nw
        Real*8  :: rhs(2)                                         ! 1 - x dimension, 2 - y dimension
    End Type Stencil

    Type(Stencil), Allocatable  ::  Qr(:)
    Real*8       , Allocatable  ::  f1(:), f2(:)
    Real                        ::  f1_n, f1_s, f2_w, f2_e, SorA  ! Sorenson Variables used for calculating f1, f2 in IMode = 3

contains
    !**************************************************************************
    ! Used to Initialize Node Coordinates. We "guess" the position of all the
    ! nodes based on the west and south boundary nodes. "a" and "b" variables
    ! can be used to stretch or compress the generated grid in direction x or y
    ! a,b > 1 result in a stretched grid. Use a=1, b=1 to bypass this option.
    !**************************************************************************
    Subroutine InitializeMesh ( a , b )
       Use System,     Only: MyGrid, Node, INode, NOP
       Implicit None
       Integer         ::   i, j
       Real*8          ::   StepX, StepY
       Real            ::   a, b
       Real*8          ::   NodeRx( MyGrid%NNodes(2) - 2 ), NodeRy( MyGrid%NNodes(1) - 2 )
       Intent(IN)      ::   a, b


       ! Define a step vector for building the initial Mesh
       StepX  =  a*ABS( Node(1)%R(1) - Node( MyGrid%TotalNodes - MyGrid%NNodes(2) + 1  )%R(1) ) / MyGrid%NNodes(1)
       StepY  =  b*ABS( Node(1)%R(2) - Node(             MyGrid%NNodes(2)              )%R(2) ) / MyGrid%NNodes(2)


       ! Assign Initial Value to Node Coordinates
       i = 1
       do INode = 2, MyGrid%NNodes(2) - 1
           NodeRx(i) = Node(INode)%R(1)
           i         = i + 1
       enddo

       INode = 1 + MyGrid%NNodes(2)
       do i = 2, MyGrid%NNodes(1) - 1
           NodeRy(i-1) = Node(INode)%R(2)
           INode       = INode + MyGrid%NNodes(2)
       enddo


       ! Generate Initial Mesh
       do j = 2, MyGrid%NNodes(1) - 1
           do i = 2, MyGrid%NNodes(2) - 1
               INode            = NOP(i,j)

               Node(INode)%R(1) = NodeRx(i-1) +  (j-1)*StepX
               Node(INode)%R(2) = NodeRy(j-1) +  (i-1)*StepY
           enddo
       enddo

    End Subroutine InitializeMesh



    !**************************************************************************
    ! Calculates Metrics of INode.
    ! IFlag: interior: 0 , South or North: 1, East or West: 2
    !**************************************************************************
    Subroutine CalculateMetrics ( INode , IFlag )
       Use System,     Only : Node, Metrics
       Implicit None
       Integer         :: IFlag
       Integer         :: INode, JNode, KNode
       Real*8          :: Metric_J, Metric_K
       Real*8          :: cov11, cov22, cov12
       Intent(IN)      :: INode, IFlag


       ! Calculate First Order Metrics
       !------------------------------
       if (IFlag /= 2 ) then
            JNode              = Node(INode)%East
            KNode              = Node(INode)%West
            Metrics(INode)%Xc  = ( Node(JNode)%R(1) - Node(KNode)%R(1) ) / 2.d0

            JNode              = Node(INode)%East
            KNode              = Node(INode)%West
            Metrics(INode)%Yc  = ( Node(JNode)%R(2) - Node(KNode)%R(2) ) / 2.d0
       endif

       if (IFlag /= 1 ) then
           JNode              = Node(INode)%North
           KNode              = Node(INode)%South
           Metrics(INode)%Xn  = ( Node(JNode)%R(1) - Node(KNode)%R(1) ) / 2.d0

           JNode              = Node(INode)%North
           KNode              = Node(INode)%South
           Metrics(INode)%Yn  = ( Node(JNode)%R(2) - Node(KNode)%R(2) ) / 2.d0
       endif

       if (IFlag == 0 ) then
            Metrics(INode)%Jac = Metrics(INode)%Xc * Metrics(INode)%Yn  -  Metrics(INode)%Xn * Metrics(INode)%Yc
            if (Metrics(INode)%Jac  < 0) then ; print*, 'Problem with negative Jacobian Value at Node:', INode ; STOP ; endif
            if (Metrics(INode)%Jac == 0) then ; print*, ' Problem with ZERO Jacobian Value at Node   :', INode ; STOP ; endif
       endif


       ! Calculate Second Order Metrics
       !--------------------------------
        if (IFlag /= 2 ) then
            JNode              = Node(INode)%East
            KNode              = Node(INode)%West
            Metrics(INode)%Xcc = ( Node(JNode)%R(1) - 2*Node(INode)%R(1) + Node(KNode)%R(1) ) / 4.d0

            JNode              = Node(INode)%East
            KNode              = Node(INode)%West
            Metrics(INode)%Ycc = ( Node(JNode)%R(2) - 2*Node(INode)%R(2) + Node(KNode)%R(2) ) / 4.d0
        endif

        if (IFlag /= 1 ) then
            JNode              = Node(INode)%North
            KNode              = Node(INode)%South
            Metrics(INode)%Xnn = ( Node(JNode)%R(1) - 2*Node(INode)%R(1) + Node(KNode)%R(1) ) / 4.d0

            JNode              = Node(INode)%North
            KNode              = Node(INode)%South
            Metrics(INode)%Xnn = ( Node(JNode)%R(2) - 2*Node(INode)%R(2) + Node(KNode)%R(2) ) / 4.d0
        endif

        if (IFlag == 0 ) then
            JNode              = Node(INode)%North
            KNode              = Node(JNode)%West
            JNode              = Node(JNode)%East
            Metric_J           = ( Node(JNode)%R(1) - Node(KNode)%R(1) ) / 2.d0
            JNode              = Node(INode)%South
            KNode              = Node(JNode)%West
            JNode              = Node(JNode)%East
            Metric_K           = ( Node(JNode)%R(1) - Node(KNode)%R(1) ) / 2.d0
            Metrics(INode)%Xcn = ( Metric_J - Metric_K )                 / 2.d0

            JNode              = Node(INode)%East
            KNode              = Node(JNode)%South
            JNode              = Node(JNode)%North
            Metric_J           = ( Node(JNode)%R(1) - Node(KNode)%R(1) ) / 2.d0
            JNode              = Node(INode)%West
            KNode              = Node(JNode)%South
            JNode              = Node(JNode)%North
            Metric_K           = ( Node(JNode)%R(1) - Node(KNode)%R(1) ) / 2.d0
            Metrics(INode)%Xnc = ( Metric_J - Metric_K )                 / 2.d0

            JNode              = Node(INode)%North
            KNode              = Node(JNode)%West
            JNode              = Node(JNode)%East
            Metric_J           = ( Node(JNode)%R(2) - Node(KNode)%R(2) ) / 2.d0
            JNode              = Node(INode)%South
            KNode              = Node(JNode)%West
            JNode              = Node(JNode)%East
            Metric_K           = ( Node(JNode)%R(2) - Node(KNode)%R(2) ) / 2.d0
            Metrics(INode)%Ycn = ( Metric_J - Metric_K )                 / 2.d0

            JNode              = Node(INode)%East
            KNode              = Node(JNode)%South
            JNode              = Node(JNode)%North
            Metric_J           = ( Node(JNode)%R(2) - Node(KNode)%R(2) ) / 2.d0
            JNode              = Node(INode)%West
            KNode              = Node(JNode)%South
            JNode              = Node(JNode)%North
            Metric_K           = ( Node(JNode)%R(2) - Node(KNode)%R(2) ) / 2.d0
            Metrics(INode)%Ync = ( Metric_J - Metric_K )

            JNode              = Node(INode)%North
            KNode              = Node(JNode)%West
            JNode              = Node(JNode)%East
            Metric_J           = ( Node(JNode)%R(1) - Node(KNode)%R(1) ) / 2.d0
            JNode              = Node(INode)%South
            KNode              = Node(JNode)%West
            JNode              = Node(JNode)%East
            Metric_K           = ( Node(JNode)%R(1) - Node(KNode)%R(1) ) / 2.d0
            Metrics(INode)%Xcn = ( Metric_J - Metric_K )                 / 2.d0

            JNode              = Node(INode)%East
            KNode              = Node(JNode)%South
            JNode              = Node(JNode)%North
            Metric_J           = ( Node(JNode)%R(1) - Node(KNode)%R(1) ) / 2.d0
            JNode              = Node(INode)%West
            KNode              = Node(JNode)%South
            JNode              = Node(JNode)%North
            Metric_K           = ( Node(JNode)%R(1) - Node(KNode)%R(1) ) / 2.d0
            Metrics(INode)%Xnc = ( Metric_J - Metric_K )                 / 2.d0

            JNode              = Node(INode)%North
            KNode              = Node(JNode)%West
            JNode              = Node(JNode)%East
            Metric_J           = ( Node(JNode)%R(2) - Node(KNode)%R(2) ) / 2.d0
            JNode              = Node(INode)%South
            KNode              = Node(JNode)%West
            JNode              = Node(JNode)%East
            Metric_K           = ( Node(JNode)%R(2) - Node(KNode)%R(2) ) / 2.d0
            Metrics(INode)%Ycn = ( Metric_J - Metric_K )                 / 2.d0

            JNode              = Node(INode)%East
            KNode              = Node(JNode)%South
            JNode              = Node(JNode)%North
            Metric_J           = ( Node(JNode)%R(2) - Node(KNode)%R(2) ) / 2.d0
            JNode              = Node(INode)%West
            KNode              = Node(JNode)%South
            JNode              = Node(JNode)%North
            Metric_K           = ( Node(JNode)%R(2) - Node(KNode)%R(2) ) / 2.d0
            Metrics(INode)%Ync = ( Metric_J - Metric_K )                 / 2.d0


            ! Calculate Covariants and overwrite them with the contravariants
            !----------------------------------------------------------------
            cov11  =  Metrics(INode)%Xc * Metrics(INode)%Xc + Metrics(INode)%Yc * Metrics(INode)%Yc
            cov22  =  Metrics(INode)%Xn * Metrics(INode)%Xn + Metrics(INode)%Yn * Metrics(INode)%Yn
            cov12  =  Metrics(INode)%Xc * Metrics(INode)%Xn + Metrics(INode)%Yc * Metrics(INode)%Yn

            ! calculate contravariants
             Metrics(INode)%g11 = + cov22 / Metrics(INode)%Jac*Metrics(INode)%Jac
             Metrics(INode)%g22 = + cov11 / Metrics(INode)%Jac*Metrics(INode)%Jac
             Metrics(INode)%g12 = - cov12 / Metrics(INode)%Jac*Metrics(INode)%Jac

       endif

    End Subroutine CalculateMetrics


    !**************************************************************************
    ! Solve a set of equations using MSIP
    !**************************************************************************
    Subroutine MSIP_Interface ( IMode , psi , Residual )
       Use System,                 Only: MyGrid
       Implicit None
       Integer                     ::  IMode
       Real*8                      ::  psi, Residual
       Intent(IN)                  ::  psi, IMode
       Intent(OUT)                 ::  Residual

       Call CreateStencil ( IMode )

       call MSIP9 ( psi , MyGrid%TotalNodes , Residual )

       Deallocate ( Qr, f1, f2 )

    End Subroutine MSIP_Interface



    !**************************************************************************
    ! Generates the Matrix to be solved based on metrics and f1, f2
    !**************************************************************************
    Subroutine CreateStencil ( IMode )
       Use System,      Only: MyGrid, Node, Metrics
       Implicit None
       Integer          ::  INode, i, IMode
       Intent(IN)       ::  IMode

       Allocate ( Qr( MyGrid%TotalNodes ) )
       Allocate ( f1( MyGrid%TotalNodes ) , f2( MyGrid%TotalNodes ) )

       do INode = 1,  MyGrid%TotalNodes
           Qr(INode)%n      =  0.d0
           Qr(INode)%s      =  0.d0
           Qr(INode)%w      =  0.d0
           Qr(INode)%e      =  0.d0
           Qr(INode)%ne     =  0.d0
           Qr(INode)%nw     =  0.d0
           Qr(INode)%se     =  0.d0
           Qr(INode)%sw     =  0.d0
           Qr(INode)%p      =  0.d0

           Qr(INode)%rhs(:) =  0.d0
       enddo

       ! Calculate f1, f2
       Call Calculate_f1_f2 ( IMode )

       do INode = 1, MyGrid%TotalNodes
           Qr(INode)%n  = + Metrics(INode)%g22 + 0.5d0 * f2(INode) * Metrics(INode)%Jac*Metrics(INode)%Jac
           Qr(INode)%s  = + Metrics(INode)%g22 - 0.5d0 * f2(INode) * Metrics(INode)%Jac*Metrics(INode)%Jac
           Qr(INode)%w  = + Metrics(INode)%g11 - 0.5d0 * f1(INode) * Metrics(INode)%Jac*Metrics(INode)%Jac
           Qr(INode)%e  = + Metrics(INode)%g11 + 0.5d0 * f1(INode) * Metrics(INode)%Jac*Metrics(INode)%Jac
           Qr(INode)%ne = + Metrics(INode)%g12 / 2.d0
           Qr(INode)%nw = - Metrics(INode)%g12 / 2.d0
           Qr(INode)%se = - Metrics(INode)%g12 / 2.d0
           Qr(INode)%sw = + Metrics(INode)%g12 / 2.d0
           Qr(INode)%p  = + 2 * (  Metrics(INode)%g11 +  Metrics(INode)%g22 )
       enddo

       ! West Boundary Conditions
       INode = 1
       do i = 1, MyGrid%NNodes(2)
           Qr(INode)%n      =  0.d0
           Qr(INode)%s      =  0.d0
           Qr(INode)%w      =  0.d0
           Qr(INode)%e      =  0.d0
           Qr(INode)%ne     =  0.d0
           Qr(INode)%nw     =  0.d0
           Qr(INode)%se     =  0.d0
           Qr(INode)%sw     =  0.d0
           Qr(INode)%p      =  1.d0

           Qr(INode)%rhs(:) =  Node(INode)%R(:)
           INode            =  INode + 1
       enddo

       ! North Boundary Conditions
       INode = 2*MyGrid%NNodes(2)
       do i = 2, MyGrid%NNodes(1)
           Qr(INode)%n      =  0.d0
           Qr(INode)%s      =  0.d0
           Qr(INode)%w      =  0.d0
           Qr(INode)%e      =  0.d0
           Qr(INode)%ne     =  0.d0
           Qr(INode)%nw     =  0.d0
           Qr(INode)%se     =  0.d0
           Qr(INode)%sw     =  0.d0
           Qr(INode)%p      =  1.d0

           Qr(INode)%rhs(:) =  Node(INode)%R(:)
           INode     = INode + MyGrid%NNodes(2)
       enddo

       ! East Boundary Conditions
       INode = MyGrid%TotalNodes - 1
       do i = 2, MyGrid%NNodes(2)
           Qr(INode)%n      =  0.d0
           Qr(INode)%s      =  0.d0
           Qr(INode)%w      =  0.d0
           Qr(INode)%e      =  0.d0
           Qr(INode)%ne     =  0.d0
           Qr(INode)%nw     =  0.d0
           Qr(INode)%se     =  0.d0
           Qr(INode)%sw     =  0.d0
           Qr(INode)%p      =  1.d0

           Qr(INode)%rhs(:) =  Node(INode)%R(:)
           INode            =  INode - 1
       enddo

      ! South Boundary Conditions
      INode = MyGrid%TotalNodes - 2*MyGrid%NNodes(2) + 1
      do i = 2, MyGrid%NNodes(1) - 1
           Qr(INode)%n      =  0.d0
           Qr(INode)%s      =  0.d0
           Qr(INode)%w      =  0.d0
           Qr(INode)%e      =  0.d0
           Qr(INode)%ne     =  0.d0
           Qr(INode)%nw     =  0.d0
           Qr(INode)%se     =  0.d0
           Qr(INode)%sw     =  0.d0
           Qr(INode)%p      =  1.d0

           Qr(INode)%rhs(:) =  Node(INode)%R(:)
           INode            = INode - MyGrid%NNodes(2)
      enddo

   End Subroutine CreateStencil

   !**************************************************************************
   ! Calculate f1, f2 based on input Mode
   !**************************************************************************
    Subroutine Calculate_f1_f2 ( IMode )
       Use System,         Only: MyGrid, Metrics, NOP
       Implicit None
       Integer             :: IMode, INode, nv, nh, i, j
       Real*8, Allocatable :: phi_n(:), phi_s(:), psi_e(:), psi_w(:)   ! Thomas-Middlecoff variables
       Real*8              :: phi, psi                                 ! Thomas-Middlecoff variables
       Real*8              :: a, b
       Intent(IN)          :: IMode

       ! Initialize f1, f2
       f1(:)   =   0.d0
       f2(:)   =   0.d0

       if ( (IMode.NE.1).AND.(IMode.NE.2).AND.(IMode.NE.3) ) then
           write(*,*)
           write(*,*) '**********************!!!!!!***********************'
           write(*,*) ' MODE ERROR: WRONG MODE # - CAN NOT FIND THIS MODE '
           write(*,*) '**********************!!!!!!***********************'
           STOP
       endif


       Select Case ( IMode )

       ! Laplace
       ! ----------------------------------
       case ( 1 )
           ! Do Nothing
           ! Use initialized f1, f2 = 0


       ! Thomas-Middlecoff
       ! ----------------------------------
       case(2)
           nv       = MyGrid%NNodes(1) - 1
           Allocate ( phi_n(2:nv) , phi_s(2:nv) )

           nh       = MyGrid%NNodes(2) - 1
           Allocate ( psi_e(2:nh) , psi_w(2:nh) )

           ! Phi North
           INode = 2*MyGrid%NNodes(2)
           do i = 2, MyGrid%NNodes(1) - 1
              call         CalculateMetrics( INode , 1 )
              a            = Metrics(INode)%Yc*Metrics(INode)%Ycc + Metrics(INode)%Xc*Metrics(INode)%Xcc
              b            = Metrics(INode)%Xc*Metrics(INode)%Xc  + Metrics(INode)%Yc*Metrics(INode)%Yc
              phi_n(i)     = - a / b
              INode        = INode + MyGrid%NNodes(2)
           enddo

           ! Phi South
           INode = MyGrid%TotalNodes - 2*MyGrid%NNodes(2) + 1
           do i = 2, MyGrid%NNodes(1) - 1
              call         CalculateMetrics( INode , 1 )
              a            = Metrics(INode)%Yc*Metrics(INode)%Ycc + Metrics(INode)%Xc*Metrics(INode)%Xcc
              b            = Metrics(INode)%Xc*Metrics(INode)%Xc  + Metrics(INode)%Yc*Metrics(INode)%Yc
              phi_s(i)     = - a / b
              INode        = INode - MyGrid%NNodes(2)
           enddo

           ! Psi West
           INode = 2
           do i = 2, MyGrid%NNodes(2) - 1
              call         CalculateMetrics( INode , 2 )
              a            = Metrics(INode)%Yn*Metrics(INode)%Ynn + Metrics(INode)%Xn*Metrics(INode)%Xnn
              b            = Metrics(INode)%Xn*Metrics(INode)%Xn  + Metrics(INode)%Yn*Metrics(INode)%Yn
              psi_w(i)     = - a / b
              INode        = INode + 1
           enddo

           ! Psi East
           INode = MyGrid%TotalNodes - 1
           do i = 2, MyGrid%NNodes(2) - 1
              call         CalculateMetrics( INode , 2 )
              a            = Metrics(INode)%Yn*Metrics(INode)%Ynn + Metrics(INode)%Xn*Metrics(INode)%Xnn
              b            = Metrics(INode)%Xn*Metrics(INode)%Xn  + Metrics(INode)%Yn*Metrics(INode)%Yn
              psi_e(i)     = - a / b
              INode        = INode - 1
           enddo

           ! Calculate f1, f2 at interior nodes based on boundaries
           do i = 2, MyGrid%NNodes(2) - 1
              do j = 2, MyGrid%NNodes(1) - 1
                   INode  =  NOP(i,j)

                   a      =  (i - 1) / ( MyGrid%NNodes(2) - 1 )
                   phi    =  phi_s(j) * ( 1 - a ) + phi_n(j) * a
                   b      =  (j - 1) / ( MyGrid%NNodes(1) - 1 )
                   psi    =  psi_w(i) * ( 1 - b ) + psi_e(i) * b

                   f1(INode) = phi * Metrics(INode)%g11
                   f2(INode) = psi * Metrics(INode)%g22
              enddo
           enddo


       ! Sorenson
       ! ------------------
       case( 3 )

           ! f1, f2 at boundaries, imposed arbitrarily
           ! Calculate f1, f2 at interior nodes
           do i = 2, MyGrid%NNodes(2) - 1
              do j = 2, MyGrid%NNodes(1) - 1
                   INode  =  NOP(i,j)

                   f1(INode) = f1_s * EXP(-SorA*float(i)) + f1_n * EXP( SorA*float(i - MyGrid%NNodes(2)) )
                   f2(INode) = f2_w * EXP(-SorA*float(j)) + f2_e * EXP( SorA*float(j - MyGrid%NNodes(1)) )
              enddo
           enddo

       End Select

    End Subroutine Calculate_f1_f2


    !*************************************************************
    ! Input Sorenson Variables - Use this subroutine to ask for
    ! SorensOn Variables before calling MSIP_Interface to solve
    !*************************************************************
    Subroutine InputSorenson
       Implicit None

       print*,   '    ***Input Sorenson Parameters***    '
       print*,   'f1 North - f1 South - f2 West - f2 East'
       read(*,*)  f1_n , f1_s , f2_w , f2_e
       print*,   'Input Sorenson "Strength" Paramater'
       read(*,*)  SorA
       print*,   '-------------------------------------'
    End Subroutine InputSorenson



    !**************************************************************************
    ! Performs MSIP factorization and solves for grid node coordinates. Updated
    ! coordinates overwrite old system coordinates
    !**************************************************************************
    Subroutine MSIP9 ( psi , NN , Residual )
       Use System,        Only: Node
       Implicit None
       Integer                 ::  NN, INode, in, is, iw, ie, ise, isw, ine, inw
       Real*8                  ::  psi, Residual
       Real*8, Dimension(0:NN) ::  am , bm , cm , dm , em , fm , gm , hm , km ! Small   letters
       Real*8                  ::  ac , bc , cc , dc , ec , fc , gc , hc , kc ! Capital letters
       Intent(IN)              ::  psi, NN
       Intent(OUT)             ::  Residual

       ! Initialization
       ! --------------
       do INode = 1, NN
           am(INode) = 0.d0
           bm(INode) = 0.d0
           cm(INode) = 0.d0
           dm(INode) = 0.d0
           em(INode) = 0.d0
           fm(INode) = 0.d0
           gm(INode) = 0.d0
           hm(INode) = 0.d0
           km(INode) = 0.d0
       enddo

       do INode = 1, NN

           ! Set indices of neighbors
           is  = Node(INode)%South
           in  = Node(INode)%North
           iw  = Node(INode)%West
           ie  = Node(INode)%East
           ine = Node(INode)%NorthEast
           inw = Node(INode)%NorthWest
           ise = Node(INode)%SouthEast
           isw = Node(INode)%SouthWest

           ! Capital Letters
           ec = - Qr(INode)%p
           hc =   Qr(INode)%e
           bc =   Qr(INode)%w
           fc =   Qr(INode)%n
           dc =   Qr(INode)%s
           kc =   Qr(INode)%ne
           cc =   Qr(INode)%nw
           gc =   Qr(INode)%se
           ac =   Qr(INode)%sw

           ! Small Letters
           am(INode)   =   ac
           bm(INode)   =   (bc - psi*cc*fm(inw) - am(INode)*fm(isw)) / (1.d0 - psi*fm(iw)*fm(inw))
           cm(INode)   =   cc - bm(INode)*fm(iw)
           dm(INode)   =   (dc - am(INode)*(2.*psi*gm(isw) + hm(isw)) - bm(INode)*gm(iw)) / (1.d0 + 2.*psi*gm(is))
           em(INode)   =   ec + am(INode)*(psi*gm(isw) - km(isw)) - bm(INode)*hm(iw) + cm(INode)*(2.*psi*fm(inw) - &
                           gm(inw) + psi*km(inw)) + dm(INode)*(2.*psi*gm(is) - fm(is))

           if ( dabs(em(INode) ).LT.1e-11)  STOP  ' Stop in Factor'

           em(INode)   =   1.d0 / em(INode)       ! attention keeps the INVERSE of EPSILON
           fm(INode)   =   (fc - bm(INode)*km(iw) - cm(INode)*(hm(inw) + 2.*psi*fm(inw) + 2.*psi*km(inw))) * em(INode)
           gm(INode)   =   (gc - dm(INode)*hm(is)) * em(INode)
           hm(INode)   =   (hc - dm(INode)*(psi*gm(is)+km(is))) * em(INode)
           km(INode)   =   kc *  em(INode)

       Enddo
       Call Solver ( NN , am , bm , cm, dm , em, fm , gm , hm , km , Residual )

    End Subroutine MSIP9

    !**************************************************************************
    ! Solves a set of NN Equations to OVERWRITE Node Coordinates
    !**************************************************************************
    Subroutine Solver ( NN , am , bm , cm, dm , em, fm , gm , hm , km , Residual )
       Use System,       Only: Node
       Implicit None
       Integer           ::  NN , INode, JNode, ic
       Type(Stencil)     ::  Sol
       Real*8            ::  Qrr(NN)  , Residual
       Real*8            ::  am(0:NN) , bm(0:NN) , cm(0:NN) , dm(0:NN) , em(0:NN) ,&
                             fm(0:NN) , gm(0:NN) , hm(0:NN) , km(0:NN) , dsol(0:NN)
       Intent(IN)        ::  NN , am , bm , cm , dm , em , fm , gm , hm , km
       Intent(OUT)       ::  Residual


       ! solve for x and y coordinate
       Residual = 0.d0
       do ic = 1, 2
           do INode = 1, NN
               JNode  =  Node(INode)%North      ;  Sol%n  =  Node(JNode)%R(ic)
               JNode  =  Node(INode)%South      ;  Sol%s  =  Node(JNode)%R(ic)
               JNode  =  Node(INode)%East       ;  Sol%e  =  Node(JNode)%R(ic)
               JNode  =  Node(INode)%West       ;  Sol%w  =  Node(JNode)%R(ic)

               JNode  =  Node(INode)%NorthEast  ;  Sol%ne =  Node(JNode)%R(ic)
               JNode  =  Node(INode)%NorthWest  ;  Sol%nw =  Node(JNode)%R(ic)
               JNode  =  Node(INode)%SouthEast  ;  Sol%se =  Node(JNode)%R(ic)
               JNode  =  Node(INode)%SouthWest  ;  Sol%sw =  Node(JNode)%R(ic)

               Sol%p      =  Node(INode)%R(ic)

               Qrr(INode) =  - ( Qr(INode)%e*Sol%e + Qr(INode)%w*Sol%w + Qr(INode)%n*Sol%n + Qr(INode)%s*Sol%s +   &
                             Qr(INode)%sw*Sol%sw + Qr(INode)%se*Sol%se + Qr(INode)%nw*Sol%nw + Qr(INode)%ne*Sol%ne &
                             - Qr(INode)%p*Sol%p + Qr(INode)%rhs(ic) )
           enddo

           call BackFront( NN , Qrr , am , bm , cm , dm , em , fm , gm , hm , km , dsol)

           do INode = 1, NN
               Node(INode)%R(ic)  =  Node(INode)%R(ic) + dsol(INode)
               Residual = dsol(INode)*dsol(INode) + Residual
           enddo
       enddo

       Residual = SQRT( Residual ) / ( 2*NN )

    End Subroutine Solver


    !********************************************************************************
    ! Executes BackFront Substitution
    !********************************************************************************
    Subroutine BackFront ( NN , Qrr , am , bm , cm , dm , em , fm , gm , hm , km , tp)
       Use System,                 Only: Node
       Implicit None
       Integer                     ::  NN, INode, in, is, iw, ie, ise, isw, ine, inw
       Type(Stencil)               ::  com
       Real*8, Dimension(NN)       ::  Qrr
       Real*8, Dimension(0:NN)     ::  am , bm , cm , dm , em , fm , gm , hm , km , Aux , tp
       Intent(IN)                  ::  NN , am , bm , cm , dm , em , fm , gm , hm , km
       Intent(OUT)                 ::  tp


       ! Initialize Aux and tp
       Aux(:)  =  0.d0
       tp(:)   =  0.d0


       ! Front Sub
       do INode = 1, NN
           ! Set indices of neighbors
           is  = Node(INode)%South
           iw  = Node(INode)%West
           inw = Node(INode)%NorthWest
           isw = Node(INode)%SouthWest

           com%w   =  0.d0
           com%sw  =  0.d0
           com%nw  =  0.d0
           com%s   =  0.d0

           if(isw.ne.0) com%sw  =  am(INode)*aux(isw)
           if(iw.ne.0)  com%w   =  bm(INode)*aux(iw)
           if(inw.ne.0) com%nw  =  cm(INode)*aux(inw)
           if(is.ne.0)  com%s   =  dm(INode)*aux(is)

           aux(INode)  =  ( Qrr(INode) - com%sw - com%w - com%nw - com%s ) * em(INode)
       enddo

       ! Back Sub
       do INode = 1, NN
           ! Set indices of neighbors
           in  = Node(INode)%North
           ie  = Node(INode)%East
           ine = Node(INode)%NorthEast
           ise = Node(INode)%SouthEast

           com%e   =  0.d0
           com%se  =  0.d0
           com%ne  =  0.d0
           com%n   =  0.d0

           if(in.ne.0)  com%n  = fm(INode)*tp(in)
           if(ise.ne.0) com%se = gm(INode)*tp(ise)
           if(ie.ne.0)  com%e  = hm(INode)*tp(ie)
           if(ine.ne.0) com%ne = km(INode)*tp(ine)

           tp(INode)  =  Aux(INode) - com%n - com%se - com%e - com%ne
       enddo
    End Subroutine BackFront


End Module GridGen
