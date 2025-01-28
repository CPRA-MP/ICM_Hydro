!> @file
!> @brief This is the subroutine to close input/output files.
!> @details This is the subroutine to close the input/output hydro files

	Subroutine CloseFiles

	use params

	implicit none

	close(204)
	close(206)
	close(207)
	close(208)
	close(209)
	close(75)
	close(96)
	close(100)
	close(103)
	close(104)
	close(1045)
	close(105)
	close(111)
	close(112)
	close(124)
	close(210)
	close(211)
	close(212)
	close(93)
	close(301)
	close(302)
	close(303)
	close(304)
	close(305)
	close(306)
	close(307)
	close(308)
	if (iWQ>0) then
		close(70)
		close(71)
		close(72)
		close(73)
		close(91)
		close(92)
		close(94)
		close(95)
		close(97)
		close(113)
		close(119)
		close(121)
		close(123)
	endif
	

	return
	end
