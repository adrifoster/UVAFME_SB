module Plot
use Tree
use Cells

implicit none

   type PlotData
      type(CellsData),   dimension(:,:),allocatable :: cells,deadcells
      type(SpeciesData),dimension(:),allocatable :: species

      real                                       :: basal_area, mn_host_DBH 
      real                                       :: host_perc_can 
      real                                       :: beetle_plot_prob 
      real                                       :: debris, infested_BA
      integer                                    :: beetle_gen
      integer                                    :: numspecies
      integer                                    :: numtrees
      integer                                    :: fire, wind
      integer                                    :: fireCount         
      integer                                    :: windCount
      integer                                    :: beetle_wind 
   end type PlotData

contains

   subroutine compute_plot_beetle_factors(self, hours_above)
	
		!computes plot-level spruce beetle probability
		class(PlotData), intent(inout) :: self
		real,            intent(in)    :: hours_above

		integer                        :: DBH_fact, BA_fact, can_fact
		integer                        :: beetle_rating
		real                           :: logit, pilog, rand_val
		real                           :: debris_prob
		real, dimension(7)             :: prob

		!calculate plot stand rating based on
		!Schmid & Frye 1976 Stand Ratings for Spruce Beetles
		
		!mn_host_DBH = mean spruce DBH above 24.5 cm (cm)
		if (self%mn_host_DBH .lt. 30.48) then
			DBH_fact = 1
		else if (self%mn_host_DBH .lt. 40.64 .and. self%mn_host_DBH .ge. 38.48) then
			DBH_fact = 2
		else
			DBH_fact = 3
		endif

		!basal_area: plot basal area (m2/ha)
		if (self%basal_area .lt. 22.95) then
			BA_fact = 1
		else if (self%basal_area .lt. 34.43 .and. self%basal_area .ge. 22.95) then
			BA_fact = 2
		else
			BA_fact = 3
		endif

	
		!host_perc_can: percent of spruce in the canopy
		if (self%host_perc_can .lt. 50.0) then
			can_fact = 1
		else if (self%host_perc_can .lt. 65.0 .and. self%host_perc_can .ge. 50.0) then
			can_fact = 2
		else
			can_fact = 3
		endif
 
		!calculate probability for univoltine beetles
		!based on Hansen et al. 2001 Can. Entomol. 133(6):827-841
		logit = -3.954 + 0.01944*hours_above
		pilog = 1.0/(1.0 + exp(-logit))
    
		!set beetle_gen based on pilog
		rand_val = urand()
		if (rand_val .lt. pilog) then
			self%beetle_gen = 2
		else
			self%beetle_gen = 1
		endif

		!calculate plot-level susceptibility
		beetle_rating = DBH_fact + BA_fact + can_fact
		self%beetle_plot_prob = 0.75*log(float(beetle_rating)) - 0.8

		!increase based on recent wind throw
		!beetle_wind: countdown from recent windthrow (starts at 10)
		if (self%beetle_wind .ge. 1 .and. self%beetle_wind .lt. 3) then
			self%beetle_plot_prob = self%beetle_plot_prob + 0.1
		else if (self%beetle_wind .ge. 3 .and. self%beetle_wind .lt. 6) then
			self%beetle_plot_prob = self%beetle_plot_prob + 0.2
		else if (self%beetle_wind .ge. 6) then
			self%beetle_plot_prob = self%beetle_plot_prob + 0.3
		end if

		if (self%beetle_plot_prob .ge. 1.0) then
			self%beetle_plot_prob = 1.0
		else if (self%beetle_plot_prob .le. 0.001) then
			self%beetle_plot_prob = 0.0
		end if

   end subroutine compute_plot_beetle_factors
