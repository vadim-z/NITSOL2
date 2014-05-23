      subroutine nitdflts(rinpt)

      implicit none  

      double precision rinpt(8)

c
c ------------------------------------------------------------------------
c
c This is nitflts v2.0, the subroutine to fill defaults for NITSOL
c parameters controlling non-linear iterations.
c
c ------------------------------------------------------------------------
c
c Explanation: 
c
c See nitsol.f
c
c ------------------------------------------------------------------------

      include 'nitdflts.h'

      rinpt(1) = DFLT_CHOICE1_EXP
      rinpt(2) = DFLT_CHOICE2_EXP
      rinpt(3) = DFLT_CHOICE2_COEF
      rinpt(4) = DFLT_ETA_CUTOFF
      rinpt(5) = DFLT_ETA_MAX
      rinpt(6) = DFLT_ETA_FIXED
      rinpt(7) = DFLT_THMIN
      rinpt(8) = DFLT_THMAX

      end
c -------------------- end of subroutine nitdflts  --------------------------


