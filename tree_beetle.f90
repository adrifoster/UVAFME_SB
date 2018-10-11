module Tree

   implicit none

   type, extends(SpeciesData) :: TreeData
     real                     :: diam_max
     real                     :: diam_bht
     real                     :: diam_canht
     real                     :: canopy_ht
     real                     :: foret_ht
     real                     :: forska_ht
     real                     :: leaf_bm   
     real                     :: biomC
     real                     :: biomN
     real                     :: beetle_prob, CK
     integer                  :: species_index
     integer                  :: mort_count
     integer                  :: beetle_counter
     integer                  :: tree_age
     integer                  :: stressor,row_num,col_num
     logical                  :: mort_marker
     logical                  :: infested
     logical                  :: nat_death,fire_death
   end type TreeData

contains

!==================================================================
! Methods
!==================================================================

  subroutine initialize_tree(self, tree_species, si)
    class(TreeData),             intent(inout) :: self
    type(SpeciesData), optional, intent(in)    :: tree_species 
    integer,           optional, intent(in)    :: si

  ! Constructor
    self%diam_bht   =0.0
    self%diam_canht =0.0
    self%canopy_ht  =std_ht
    self%foret_ht   =1.0
    self%forska_ht  =1.0
    self%biomC      =0.0
    self%biomN      =0.0
    self%beetle_prob =0.0
    self%CK = 0.0
    self%mort_count =0
    self%row_num =0
    self%col_num =0
    self%beetle_counter =0
    self%tree_age   =0
    self%stressor   =0
    self%mort_marker=.false.
    self%infested = .false.
    self%nat_death  =.false.
    self%fire_death =.false.

    if ( present(tree_species) ) then
       self%SpeciesData=tree_species
    endif

    if ( present(si) ) then
       self%species_index=si
    else
       self%species_index =0
    endif

  end subroutine initialize_tree


  subroutine copy_tree(self,tree)
    class(TreeData), intent(inout):: self
    class(TreeData), intent(in)   :: tree

      self%diam_max      = tree%diam_max
      self%diam_bht      = tree%diam_bht
      self%diam_canht    = tree%diam_canht
      self%canopy_ht     = tree%canopy_ht
      self%foret_ht      = tree%foret_ht
      self%forska_ht     = tree%forska_ht
      self%leaf_bm       = tree%leaf_bm
      self%biomC         = tree%biomC
      self%biomN         = tree%biomN
      self%row_num       = tree%row_num
      self%col_num       = tree%col_num
      self%beetle_prob   = tree%beetle_prob
      self%CK            = tree%CK
      self%species_index = tree%species_index
      self%mort_count    = tree%mort_count
      self%beetle_counter = tree%beetle_counter
      self%tree_age      = tree%tree_age
      self%mort_marker   = tree%mort_marker
      self%infested      = tree%infested
      self%nat_death     = tree%nat_death
      self%fire_death    = tree%fire_death
      self%stressor      = tree%stressor
      self%SpeciesData   = tree%SpeciesData

  end subroutine copy_tree

  subroutine wind_survival(tree, wcat, w_survive)  
	!calculates mortality of a tree from a windthrow event
	
	class(TreeData), intent(in)  :: tree
	real,            intent(in)  :: wcat !intensity of wind (1-10)
	logical,         intent(out) :: w_survive
	real                         :: pwind, logit, rand_val
	!diam_bht: tree DBH (cm)

	!from Rich et al 2007 J. of Ecol. 95:1261-1273
	if (wcat .ge. 0.1) then
     
		logit = 0.75*log(tree%diam_bht) 
     
		pwind = 1.0/(1.0 + exp(-logit))
     
		rand_val = urand()
     
		if (rand_val .lt. pwind) then
			w_survive = .false.
		else
			w_survive = .true.
		endif
	
	else 
		w_survive = .true.
	endif


  end subroutine wind_survival

  subroutine beetle_infested(tree, beetle_on, inf_ba)
	!determines whether tree is infested with spruce beetles
	
	class(TreeData), intent(inout) :: tree
    real,            intent(in)    :: inf_ba !infested basal area on plot (m2/ha)
    logical,         intent(in)    :: beetle_on !logical, is submodel turned on
    real                           :: rand_val, tree_min

    rand_val = urand()
    
    if (beetle_on) then

		if (inf_ba .ge. 10.0) then
			!minimum tree size that can be infested is 10 cm
             tree_min = 10.0
		else
			!minimum tree size that can be infested is 30 cm
             tree_min = 30.0
		endif
	
		if ((tree%beetle_host) .and. (tree%diam_bht .gt. tree_min)) then
        
			if (rand_val .lt. tree%beetle_prob) then
				tree%infested = .true.
			else
				tree%infested = .false. !tree is not infested
        
			endif
      
		else
			tree%infested = .false. !tree is not a host or is too small
		endif
    
	else
    
		tree%infested = .false. !beetle submodel not turned on	
	endif

  end subroutine beetle_infested

  subroutine beetle_survival(tree, beetle_on, b_survive)
	!determines if tree has died from spruce beetle infestation
    
    class(TreeData), intent(in)  :: tree
    logical,         intent(in)  :: beetle_on !is submodel turned on
    logical,         intent(out) :: b_survive

	if (beetle_on) then
		if (tree%beetle_host) then
			if (tree%infested) then   
				if (tree%beetle_counter .ge. 5) then
					!trees remain standing dead after infestation for
					!5 years
					b_survive = .false.
				
				else
					b_survive = .true. !hasn't been 5 years yet
				endif
			else
				b_survive = .true. !tree is not infested
			endif
		else
			b_survive = .true. !tree is not a host
		endif
	else
		b_survive = .true. !beetle submodel not turned on
    endif
    
  end subroutine beetle_survival
 
     
  subroutine beetle_prob_calc(tree, site_prob, beetle_on, beetle_gen,  &
						debris, inf_ba)
  
	  !calculates probability of tree being infested by spruce beetles
      class(TreeData), intent(inout) :: tree
      real,            intent(in)    :: site_prob, debris, inf_ba
      integer,         intent(in)    :: beetle_gen
      logical,         intent(in)    :: beetle_on
      real                           :: DBH_prob, stress_prob, age_prob
      real                           :: CK_prob, debris_prob, tree_min
      real                           :: plot_rand, init_prob
      real                           :: beetle_prob, gen_fact
      
      !inf_ba: current year's plot-level infested basal area (m2/ha)
      !site_prob: probability of infestation from site charactersitics
      !beetle_on: is submodel turned on
      !beetle_gen: semivoltine or univoltine beetles on plot
		!1: semivoltine beetles
		!2: univoltine beetles
	  !debris: amount of coarse woody debris on plot (tonnes C/ha)
      
      if (beetle_on) then
         
         		if (tree%beetle_host) then

          			if (inf_ba .ge. 10.0) then
             			tree_min = 10.0 !minimum infestable DBH (cm)
          			else
            			tree_min = 30.0 !minimum infestable DBH (cm)
          			endif

         			if (tree%diam_bht .gt. tree_min) then
						
						!calculate stress-induced susceptibility to infestation
						!based on how many years of prolonged low-diameter
						!increment growth (mort_count) of tree
             			if (tree%mort_count .eq. 0) then
               				stress_prob = 0.0
             			elseif (tree%mort_count .eq. 1) then
               				stress_prob = 0.1
             			elseif (tree%mort_count .eq. 2) then
               				stress_prob = 0.2
             			elseif (tree%mort_count .eq. 3) then
               				stress_prob = 0.3
             			endif
      
             			!calculate susceptibility from DBH (cm)	
             			DBH_prob = 0.010521*tree%diam_bht
             			if (DBH_prob .le. 0.0001) then
               				DBH_prob = 0.0
             			else if (DBH_prob .ge. 1.0) then
               				DBH_prob = 1.0
             			end if
           
						!calculate susceptibility from damage from fires
						!CK = scrown scorch (0-100%) from 
						!Keane et al. 2011 USDA Forest Service Report
						!RMRS-GTR-55:145
						
             			CK_prob = tree%CK/100 
             			if (CK_prob .le. 0.001) then
                			CK_prob = 0.0
             			elseif(CK_prob .ge. 1.0) then
                			CK_prob = 1.0
             			endif

						!calculate susceptibility from amount of CWD on
						!plot
             			debris_prob = (debris*20.0)/300.0
             			if (debris_prob .ge. 1.0) then
               				debris_prob = 1.0
             			elseif (debris_prob .le. 0.0001) then
               				debris_prob = 0.0
             			endif


						!calculate tree-level susceptibility
						!adapted from Seidl et al. 2007 Ecol. Model. 206(3-4):383-399
             			init_prob = min((0.5*site_prob + 0.2*stress_prob +       &
                           	0.25*DBH_prob + 0.1*CK_prob + 0.4*debris_prob), 1.0)             

						!set gen_fact based on univoltine or semivoltine
						!beetles on plot
             			if (beetle_gen .eq. 2) then
							!univoltine beetles
                 			gen_fact = 1.8
             			else
							!semivoltine beetles
                 			gen_fact = 0.5
            	 		endif
           
						!calculate final beetle infestation probability
            	 		beetle_prob = 1 - exp(-2.0*init_prob**1.3)**gen_fact

             			if (beetle_prob .gt. 1.0) then
                			tree%beetle_prob = 1.0
             			else if (beetle_prob .lt. 0.0001) then
                			tree%beetle_prob = 0.0
             			else 
                			tree%beetle_prob = beetle_prob   
             			endif

         		 	else 
             			tree%beetle_prob = 0.0 !tree is too small to be infested
          			endif

       			 else
           			tree%beetle_prob = 0.0  !tree is not a host
        		 endif

				else
					tree%beetle_prob = 0.0 !plot not eligible
				end if
      	else
        		tree%beetle_prob = 0.0  !beetle submodel turned off
      	endif

  end subroutine beetle_prob_calc


end module Tree
