Module System
   Implicit None

   Type Mesh
      Integer   :: TotalNodes                                 ! Total # of Nodes
      Integer   :: NNodes(2)                                  ! # of Cells in each direction
   End Type Mesh

   Type Nodes
      Integer  :: IJ(2)                                       ! (I, J) Indices
      Integer  :: North, South, West, East                    ! Neighbors
      Integer  :: NorthEast, NorthWest, SouthEast, SouthWest  ! Neighbors
      Real*8   :: R(2)                                        ! Coordinates x, y
   End Type Nodes

   Type Metric
      Real*8   :: Xc, Yc, Xn, Yn                              ! First Order Metrics
      Real*8   :: Xcc, Ycc, Xnn, Ynn, Xcn, Xnc, Ycn, Ync      ! Second Order Metrics
      Real*8   :: g11, g12, g22                               ! Contravariants
      Real*8   :: Jac                                         ! Jacobian
   End Type Metric

   Integer                    :: INode                        ! Global Node #
   Integer     , Allocatable  :: NOP(:,:)                     ! NOP array
   Type(Mesh)                 :: MyGrid
   Type(Nodes) , Allocatable  :: Node(:)                      ! Node(INode) returns info on a certain node
   Type(Metric), Allocatable  :: Metrics(:)                   ! Metrics(INode) stores metric values

contains



   !***********************************************************************
   ! Reads boundary nodes from the input file and allocates module arrays
   ! accordingly. Also calls internally subroutine string to generate node
   ! numbering.
   !***********************************************************************
   Subroutine ReadInput
      Implicit None
      Integer    ::  i

      ! Open Input File and read number of nodes in each direction (x,y)
      open( 4 , File='input' , status='UNKNOWN' )
      read(4,*) MyGrid%NNodes(1), MyGrid%NNodes(2)


      ! Calculate Total # of Nodes and allocate system variables
      MyGrid%TotalNodes  =  PRODUCT(MyGrid%NNodes(:))

      Allocate ( Node(0:MyGrid%TotalNodes) , Metrics(MyGrid%TotalNodes) )
      Allocate ( NOP( MyGrid%NNodes(2)     ,         MyGrid%NNodes(1) ) )

      ! Call Subroutine String
      Call String


      !Initialize Node Coordinates
      do INode = 0, MyGrid%TotalNodes
          Node(INode)%R(:) = 0.d0
      enddo


      ! Read West Boundary Node Coordinates Clockwise
      INode = 1
      do i = 1, MyGrid%NNodes(2)
         read(4,*) Node(INode)%R(:)
         INode            = INode + 1
      enddo

      ! Read North Boundary Node Coordinates Clockwise
      INode = 2*MyGrid%NNodes(2)
      do i = 2, MyGrid%NNodes(1)
         read(4,*) Node(INode)%R(:)
         INode     = INode + MyGrid%NNodes(2)
      enddo

      ! Read East Boundary Node Coordinates Clockwise
      INode = MyGrid%TotalNodes - 1
      do i = 2, MyGrid%NNodes(2)
         read(4,*) Node(INode)%R(:)
         INode     = INode - 1
      enddo

      ! Read South Boundary Node Coordinates Clockwise
      INode = MyGrid%TotalNodes - 2*MyGrid%NNodes(2) + 1
      do i = 2, MyGrid%NNodes(1) - 1
         read(4,*) Node(INode)%R(:)
         INode     = INode - MyGrid%NNodes(2)
      enddo

      close(4)
   End Subroutine ReadInput

   !***********************************************************************
   ! Subroutine ReadBoundaryDisplacement reads an input displacement file
   ! and updates Boundary Node coordinates. The input file should be created
   ! in a clockwise format, starting at INode = 1
   !***********************************************************************
   Subroutine ReadBoundaryDisplacement
      Implicit None
      Integer     ::  i
      Real*8      ::  Dis( 2 , MyGrid%TotalNodes )     ! Node displacement

      ! Initialize Dis - Internal Nodes will remain with ZERO displacement
      Dis(:,:) = 0.d0

      ! Open Input File
      open( 11 , File='BoundaryDis' , status='UNKNOWN' )


      ! Read West Boundary Node Displacement Clockwise
      INode = 1
      do i = 1, MyGrid%NNodes(2)
         read(11,*) Dis( : , INode )
         INode      = INode + 1
      enddo

      ! Read North Boundary Node Coordinates Clockwise
      INode = 2*MyGrid%NNodes(2)
      do i = 2, MyGrid%NNodes(1)
         read(11,*) Dis( : , INode )
         INode      = INode + MyGrid%NNodes(2)
      enddo

      ! Read East Boundary Node Coordinates Clockwise
      INode = MyGrid%TotalNodes - 1
      do i = 2, MyGrid%NNodes(2)
         read(11,*) Dis( : , INode )
         INode      = INode - 1
      enddo

      ! Read South Boundary Node Coordinates Clockwise
      INode = MyGrid%TotalNodes - 2*MyGrid%NNodes(2) + 1
      do i = 2, MyGrid%NNodes(1) - 1
         read(11,*) Dis( : , INode )
         INode      = INode - MyGrid%NNodes(2)
      enddo

      close(11)


      ! Update Node coordinates
      do INode = 1, MyGrid%TotalNodes
         Node(INode)%R(:) = Node(INode)%R(:) + Dis( : , INode )
      enddo


   End Subroutine ReadBoundaryDisplacement



   !***********************************************************************
   ! Subroutine String generates node numbering, NOP matrix and assigns
   ! neighbor nodes(north, west, east, south, northeast, etc) to each node.
   ! i -> rows   , from bottom to top
   ! j -> columns, from left to right
   !***********************************************************************
   Subroutine String
      Implicit None
      Integer    :: i, j, nx, ny


      ! Initialize Neighbors and (I,J) coordinates
      do INode = 0, MyGrid%TotalNodes
         Node(INode)%North      = 0
         Node(INode)%South      = 0
         Node(INode)%West       = 0
         Node(INode)%East       = 0
         Node(INode)%NorthEast  = 0
         Node(INode)%NorthWest  = 0
         Node(INode)%SouthEast  = 0
         Node(INode)%SouthWest  = 0

         Node(INode)%IJ(:)      = 0
      enddo

      nx = MyGrid%NNodes(1)
      ny = MyGrid%NNodes(2)

      ! Assign (I,J) Coordinates and Generate NOP Matrix
      INode = 0
      do j = 1, nx
         do i = 1, ny
            INode             = INode + 1

            Node(INode)%IJ(:) = (/ i , j /)
            NOP(i,j)          = INode
         enddo
      enddo

      ! Assign Neigbors to each Node
      do INode = 1, MyGrid%TotalNodes
         i  = Node(INode)%IJ(1)
         j  = Node(INode)%IJ(2)

         if (j.lt.nx)   Node(INode)%East   = NOP(   i  , j+1 )
         if (j.gt.1)    Node(INode)%West   = NOP(   i  , j-1 )
         if (i.lt.ny)   Node(INode)%North  = NOP( i+1  ,   j )
         if (i.gt.1)    Node(INode)%South  = NOP( i-1  ,   j )
         if ((j.lt.nx).and.(i.gt.1) )  Node(INode)%SouthEast = NOP( i-1 , j+1 )
         if ((j.gt.1) .and.(i.gt.1) )  Node(INode)%SouthWest = NOP( i-1 , j-1 )
         if ((j.gt.1) .and.(i.lt.ny))  Node(INode)%NorthWest = NOP( i+1 , j-1 )
         if ((j.lt.nx).and.(i.lt.ny))  Node(INode)%NorthEast = NOP( i+1 , j+1 )

      enddo
   End Subroutine String


   !***********************************************************************
   ! Case IFlag = 0: Deallocates arrays, Case IFlag = 1: everything ZERO
   !***********************************************************************
   Subroutine DelSystem ( IFlag )
      Implicit None
      Integer    ::  IFlag
      Intent(IN) ::  IFlag

      if ( ((IFlag /= 0).AND.(IFlag /= 1)) ) then ; print*, 'DelSystem ERROR: IFlag must be 0 or 1' ; STOP ; endif

      Select Case ( IFlag )

      Case ( 0 )
           deallocate ( NOP ) ; deallocate ( Node ) ; deallocate ( Metrics )

      Case ( 1 )
           ! Nullify NOP Matrix
           NOP(:,:) = 0

           do INode = 1, MyGrid%TotalNodes
               ! Nullify Nodes
               Node(INode)%IJ(:)     = 0
               Node(INode)%R(:)      = 0.d0
               Node(INode)%North     = 0 ; Node(INode)%South     = 0 ; Node(INode)%West      = 0 ; Node(INode)%East      = 0
               Node(INode)%NorthWest = 0 ; Node(INode)%SouthWest = 0 ; Node(INode)%SouthEast = 0 ; Node(INode)%NorthEast = 0

               ! Nullify Metrics
               Metrics(INode)%Xc  = 0.d0 ; Metrics(INode)%Yc  = 0.d0 ; Metrics(INode)%Xn  = 0.d0 ; Metrics(INode)%Yn  = 0.d0
               Metrics(INode)%Xcc = 0.d0 ; Metrics(INode)%Ycc = 0.d0 ; Metrics(INode)%Xnn = 0.d0 ; Metrics(INode)%Ynn = 0.d0
               Metrics(INode)%Xcn = 0.d0 ; Metrics(INode)%Ycn = 0.d0 ; Metrics(INode)%Xcn = 0.d0 ; Metrics(INode)%Ycn = 0.d0
               Metrics(INode)%g11 = 0.d0 ; Metrics(INode)%g12 = 0.d0 ; Metrics(INode)%g22 = 0.d0
               Metrics(INode)%Jac = 0.d0
           Enddo

      End Select

      ! Nullify MyGrid
      MyGrid%TotalNodes  = 0
      MyGrid%NNodes(:)   = 0

   End Subroutine DelSystem

   !***********************************************************
   ! Prints System Node Coordinates in MyGrid.txt or DisGrid.txt
   !************************************************************
   Subroutine PrintSystem ( IPrint )
      Implicit None
      Integer   ::  INode, IPrint
      Intent(IN)::  IPrint

      if ( (IPrint.NE.23).AND.(IPrint.NE.24) ) STOP ' !!! PRINT ERROR !!!'

      open( IPrint , File='MyGrid.txt'  , status='UNKNOWN' )
      open( IPrint , File='DisGrid.txt' , status='UNKNOWN' )

      do INode = 1, MyGrid%TotalNodes
          write(IPrint,*) INode, Node(INode)%R(:)
      enddo

      close ( IPrint )

   End Subroutine PrintSystem

End Module System
