module Main_Model

  use Soil
  use Site
  use Species
  use tree_beetle
  use hours_above
  use Cells

  implicit none

contains

  subroutine Beetle_Model(site, year)
	!calculates relevant beetle metrics on the plot
	!calculates spruce-beetle infestation probability and status

    type(SiteData),                intent(inout) :: site
    integer,                       intent(in)    :: year

    real                                         :: beetle_DBH, beetle_LAI
    real                                         :: other_LAI
    real                                         :: basal_ar, inf_BA
    integer                                      :: beetle_ntree
    integer                                      :: ntrees
    integer                                      :: is, ip, k, r, c
	
	do ip = 1, site%numplots
		
		!get number of trees
		ntrees = count(site%plots(ip)%cells(:,:)%filled)
		
		!set values to 0 for each plot
		beetle_DBH = 0.0   !average DBH of adult spruce (cm)
		beetle_LAI = 0.0   !leaf biomass of adult spruce 
		other_LAI = 0.0    !leaf biomass of adult non-host
		basal_ar = 0.0     !spruce basal area (m2)
		beetle_ntree = 0.0 !number of adult hosts
		inf_BA = 0.0       !infested basal area (m2)
	
		!loop through trees and get values for plot-susceptibility factors
		if (ntrees .gt. 0) then
			
			do r = 1, maxcells
				do c = 1, maxcells
					if (site%plots(ip)%cells(r,c)%filled) then
					if (site%beetle_on) then
					
						if (site%plots(ip)%cells(r,c)%tree%infested) then
							inf_BA = inf_BA + (0.000025*pi*            &
								site%plots(ip)%cells(r,c)%tree%diam_bht**2)
						endif
						
						if (site%plots(ip)%cells(r,c)%tree%beetle_host) then
							
							if (site%plots(ip)%cells(r,c)%tree%diam_bht .gt. 25.4) then
					
								beetle_DBH = beetle_DBH +              &
									site%plots(ip)%cells(r,c)%tree%diam_bht
								beetle_ntree = beetle_ntree + 1
								beetle_LAI = beetle_LAI +              &
									lai_biomass_c(site%plots(ip)%cells(r,c)%tree)
							endif
							
							basal_ar = basal_ar +(0.000025*pi*         &
								site%plots(ip)%cells(r,c)%tree%diam_bht**2)
							
						else
							
							if (site%plots(ip)%cells(r,c)%tree%diam_bht .gt. 25.4) then
								other_LAI = other_LAI +                &
								  lai_biomass_c(site%plots(ip)%cells(r,c)%tree)
							endif

						endif
			 
					end if
					end if
				end do
			end do
			
			! aggregate plot-level beetle factors and calculate overall plot factor
			if (site%beetle_on) then
				
				!basal area of spruce (m2/ha)
				site%plots(ip)%basal_area = basal_ar*hec_to_m2/plotsize
				site%plots(ip)%infested_BA= inf_BA*hec_to_m2/plotsize
				
				if (beetle_ntree .lt. 0.0001) then
					site%plots(ip)%mn_host_DBH = 0.0
				else
					site%plots(ip)%mn_host_DBH = beetle_DBH/float(beetle_ntree)
				endif

				if ((beetle_LAI + other_LAI) .lt. 0.0001) then
					site%plots(ip)%host_perc_can = 0.0
				else
					site%plots(ip)%host_perc_can = beetle_LAI/(beetle_LAI + other_LAI)*100
				endif
				 
			   
				call compute_plot_beetle_factors(site%plots(ip),site%hours_above)

			
				!calculate tree-level probability and infestation status
				do r = 1, maxcells 
					do c = 1, maxcells
						if (site%plots(ip)%cells(r,c)%filled) then
							
							call update_tree(site%plots(ip)%cells(r,c)%tree,site%species(k))

							if (site%plots(ip)%cells(r,c)%tree%infested .eq. .false.) then

								call beetle_prob_calc(site%plots(ip)%cells(r,c)%tree,  &
									site%plots(ip)%beetle_plot_prob,site%beetle_on,    &
									 site%plots(ip)%beetle_gen,site%plots(ip)%debris,  &
									 site%plots(ip)%infested_BA)
		 
								call beetle_infested(site%plots(ip)%cells(r,c)%tree,site%beetle_on, &
										  site%plots(ip)%infested_BA)

							end if
						end if
					end do
				end do 

				!increase probability based on proximity to other trees
				do r = 1, maxcells
					do c = 1, maxcells
						if (site%plots(ip)%cells(r,c)%filled) then

							call update_tree(site%plots(ip)%cells(r,c)%tree,  &
									site%species(k))
				
							if (site%plots(ip)%cells(r,c)%tree%infested) then


							!cardinal directions get affected more

							 if (site%plots(ip)%cells(modulo(r+1-1,maxcells)+1,c)%filled) then
							   site%plots(ip)%cells(modulo(r+1-1,maxcells)+1,          &
												   c)%tree%beetle_prob =               &
							   site%plots(ip)%cells(modulo(r+1-1,maxcells)+1,          &
												   c)%tree%beetle_prob + 0.2
							 endif


							 if (site%plots(ip)%cells(modulo(r-1-1,maxcells)+1,c)%filled) then
							   site%plots(ip)%cells(modulo(r-1-1,maxcells)+1,          &
												   c)%tree%beetle_prob =               &
							   site%plots(ip)%cells(modulo(r-1-1,maxcells)+1,          &
												 c)%tree%beetle_prob + 0.2
							endif

							 if (site%plots(ip)%cells(r,modulo(c+1-1,maxcells)+1)%filled) then
							   site%plots(ip)%cells(r,                                 &
								 modulo(c+1-1,maxcells)+1)%tree%beetle_prob =          &
							   site%plots(ip)%cells(r,                                 &
								 modulo(c+1-1,maxcells)+1)%tree%beetle_prob + 0.2
							 endif

							 if (site%plots(ip)%cells(r,modulo(c-1-1,maxcells)+1)%filled) then
							   site%plots(ip)%cells(r,                                 &
								 modulo(c-1-1,maxcells)+1)%tree%beetle_prob =          &
							   site%plots(ip)%cells(r,                                 &
								  modulo(c-1-1,maxcells)+1)%tree%beetle_prob + 0.2
							 endif


							 if (site%plots(ip)%cells(modulo(r+1-1,maxcells)+1,         &
							   modulo(c+1-1,maxcells)+1)%filled) then

								site%plots(ip)%cells(modulo(r+1-1,maxcells)+1,          &
									  modulo(c+1-1,maxcells)+1)%tree%beetle_prob =      &
								site%plots(ip)%cells(modulo(r+1-1,maxcells)+1,          &
									modulo(c+1-1,maxcells)+1)%tree%beetle_prob + 0.1
							 endif

							 if (site%plots(ip)%cells(modulo(r+1-1,maxcells)+1,         &
							   modulo(c-1-1,maxcells)+1)%filled) then

							  site%plots(ip)%cells(modulo(r+1-1,maxcells)+1,            &
										 modulo(c-1-1,maxcells)+1)%tree%beetle_prob =   &
							  site%plots(ip)%cells(modulo(r+1-1,maxcells)+1,            &
									modulo(c-1-1,maxcells)+1)%tree%beetle_prob + 0.1
							 endif

							 if (site%plots(ip)%cells(modulo(r-1-1,maxcells)+1,         &
							   modulo(c+1-1,maxcells)+1)%filled) then

							  site%plots(ip)%cells(modulo(r-1-1,maxcells)+1,            &
										 modulo(c+1-1,maxcells)+1)%tree%beetle_prob =   &
							  site%plots(ip)%cells(modulo(r-1-1,maxcells)+1,            &
									modulo(c+1-1,maxcells)+1)%tree%beetle_prob + 0.1
							endif

						   if (site%plots(ip)%cells(modulo(r-1-1,maxcells)+1,          &
							  modulo(c-1-1,maxcells)+1)%filled) then

							 site%plots(ip)%cells(modulo(r-1-1,maxcells)+1,             &
										 modulo(c-1-1,maxcells)+1)%tree%beetle_prob =   &
							 site%plots(ip)%cells(modulo(r-1-1,maxcells)+1,             &
									modulo(c-1-1,maxcells)+1)%tree%beetle_prob + 0.1
						   endif

						   endif ! end if infested
						endif !end if filled
					end do 
				end do !end fifth cells loop

			   ! recalculate infested trees based on proximity increases
				do r = 1,maxcells
					do c = 1,maxcells
						if(site%plots(ip)%cells(r,c)%filled) then

							call update_tree(site%plots(ip)%cells(r,c)%tree, &
									   site%species(k))

							if (site%plots(ip)%cells(r,c)%tree%infested .eq. .false.) then

								call beetle_infested(site%plots(ip)%cells(r,c)%tree,     &
										site%beetle_on,site%plots(ip)%infested_BA)

							endif
					
						endif ! end if filled
					enddo
				enddo !end cells loop

			end if !end if beetle_on
		end if      

    end do

  end subroutine Beetle_Model


  subroutine Beetle_Mortality(site,year)
	!kills trees that die from beetle mortality

    type(SiteData),        intent(inout) :: site
    integer,               intent(in)    :: year

    real                                 :: leaf_b
    real                                 :: leaf_bm, bmc
    integer                              :: num_species
    integer                              :: jt, k          
    integer                              :: ip, is, ht, r, c

    num_species=size(site%species)

    leaf_b = 1.0 + con_leaf_ratio

    do ip=1,site%numplots
		
		site%plots(ip)%deadcells(:,:)%filled = .false.

		site%plots(ip)%d_type = 0.0

		site%plots(ip)%numtrees = count(site%plots(ip)%cells(:,:)%filled)

		if (site%plots(ip)%numtrees > 0) then

       

			do r = 1,maxcells
				do c = 1,maxcells
					
					if (site%plots(ip)%cells(r,c)%filled) then
					
						k=site%plots(ip)%cells(r,c)%tree%species_index
						call update_tree(site%plots(ip)%cells(r,c)%tree,site%species(k))

						!spruce trees lose their leaves after 2 years
						if (site%plots(ip)%cells(r,c)%tree%beetle_counter .gt. 2) then
							leaf_bm = 0.0
						else
							call leaf_biomass_c(site%plots(ip)%cells(r,c)%tree)
							leaf_bm=site%plots(ip)%cells(r,c)%tree%leaf_bm
						endif

						if (site%plots(ip)%cells(r,c)%tree%infested .eq. .false.) then
							call growth_survival(site%plots(ip)%cells(r,c)%tree, growth_survive)
 
							call age_survival(site%plots(ip)%cells(r,c)%tree, age_survive)

						end if
							call beetle_survival(site%plots(ip)%cells(r,c)%tree,     &
								site%beetle_on,beetle_survive)

						if (growth_survive .and. age_survive .and. beetle_survive) then 

							site%plots(ip)%cells(r,c)%filled = .true.
							site%plots(ip)%cells(r,c)%tree%tree_age=   &
								site%plots(ip)%cells(r,c)%tree%tree_age + 1

							if (site%plots(ip)%cells(r,c)%tree%infested) then
									site%plots(ip)%cells(r,c)%tree%beetle_counter =        &
									site%plots(ip)%cells(r,c)%tree%beetle_counter + 1
							endif


							if (site%plots(ip)%cells(r,c)%tree%beetle_counter .eq. 2) then
							 site%soil%C_into_A0=site%soil%C_into_A0 + leaf_bm
							 site%soil%N_into_A0=site%soil%N_into_A0 +             &
											  leaf_bm/con_leaf_c_n
							

							else

								if (site%species(k)%conifer) then
									site%soil%C_into_A0=leaf_bm*(leaf_b-1.0)+              &
														site%soil%C_into_A0
									site%soil%N_into_A0=site%soil%N_into_A0+               &
											leaf_bm*(leaf_b-1.0)/con_leaf_c_n

						   
								else
									site%soil%C_into_A0=leaf_bm+site%soil%C_into_A0
									site%soil%N_into_A0=site%soil%N_into_A0+               &
															   leaf_bm/dec_leaf_c_n

								end if
							end if

						else
                   
							!tree dies
							site%plots(ip)%cells(r,c)%filled = .false.
							site%plots(ip)%deadcells(r,c)%filled=.true.

							call copy_tree(site%plots(ip)%deadcells(r,c)%tree, &
										site%plots(ip)%cells(r,c)%tree)

							if (growth_survive .eq. .false.) then
								site%plots(ip)%d_type(site%plots(ip)%cells(r,c)%tree%stressor) = &
								site%plots(ip)%d_type(site%plots(ip)%cells(r,c)%tree%stressor) + &
								site%plots(ip)%cells(r,c)%tree%biomC + leaf_bm

								if (site%plots(ip)%cells(r,c)%tree%beetle_host .and.      &
									site%plots(ip)%cells(r,c)%tree%diam_bht .gt. 30.0) then
									site%plots(ip)%debris = site%plots(ip)%debris +        &
									site%plots(ip)%cells(r,c)%tree%biomC
								end if


							else if (beetle_survive .eq. .false.) then
								site%plots(ip)%d_type(8) = site%plots(ip)%d_type(8) + &
									site%plots(ip)%cells(r,c)%tree%biomC + leaf_bm
								site%plots(ip)%deadcells(r,c)%tree%stressor = 8
							end if
                   
							bmc=site%plots(ip)%cells(r,c)%tree%biomC

						   if (site%species(k)%conifer) then 
							  site%soil%C_into_A0=site%soil%C_into_A0+bmc+             &
																	leaf_bm*leaf_b
							  site%soil%N_into_A0=site%soil%N_into_A0+bmc/stem_c_n +   &
														leaf_bm/con_leaf_c_n*leaf_b
		 

						   else
							  site%soil%C_into_A0=site%soil%C_into_A0+bmc+leaf_bm
							  site%soil%N_into_A0=site%soil%N_into_A0+bmc/stem_c_n+    &
															   leaf_bm/dec_leaf_c_n

						end if         
					end if
				end do 
			end do            
 
         endif               
          
          site%plots(ip)%numtrees=count(site%plots(ip)%cells(:,:)%filled)


	end do ! end plot loop


    return

  end subroutine Beetle_Mortality

end module Main_Model
