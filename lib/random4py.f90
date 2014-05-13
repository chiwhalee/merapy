module mRandom
   integer JSEED
   integer IFRST
   integer NEXTN
   integer COUNT_NUM
contains
   SUBROUTINE cSRAND(ISEED)
!!$
!!$  This subroutine sets the integer seed to be used with the
!!$  companion RAND function to the value of ISEED.  A flag is 
!!$  set to indicate that the sequence of pseudo-random numbers 
!!$  for the specified seed should start from the beginning.
!!$
!!$
      JSEED = ISEED
      NEXTN = JSEED
      IFRST = 0
      COUNT_NUM = 0
   END SUBROUTINE CSRAND
   
   REAL*8 FUNCTION cRAND()
      PARAMETER (MPLIER=16807,MODLUS=2147483647,&
           MOBYMP=127773,MOMDMP=2836)
      INTEGER HVLUE, LVLUE, TESTV!, NEXTN

      COUNT_NUM = COUNT_NUM + 1
      IF (IFRST .EQ. 0) THEN
        NEXTN = JSEED
        IFRST = 1
     ENDIF
     HVLUE = NEXTN / MOBYMP
     LVLUE = MOD(NEXTN, MOBYMP)
     TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
     IF (TESTV .GT. 0) THEN
        NEXTN = TESTV
     ELSE
        NEXTN = TESTV + MODLUS
     ENDIF
     cRAND = dble(NEXTN)/dble(MODLUS)
     RETURN
  END FUNCTION CRAND
  
end module mRandom
