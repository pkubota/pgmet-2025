MODULE Class_Module_TimeManager
 USE Constants, Only: InitClassModuleConstants, r8,r4,i4,nfprt
 IMPLICIT NONE
 PRIVATE
 INTEGER         , PUBLIC           :: idatei(4)
 INTEGER         , PUBLIC           :: idatec(4)
 INTEGER         , PUBLIC           :: idatef(4)
 INTEGER         , PUBLIC           :: DHFCT
 INTEGER         , PUBLIC           :: cth0
 REAL(KIND=r8)   , PUBLIC           :: dct
 REAL(KIND=r8)   , PUBLIC           :: dt_step
 INTEGER         , PUBLIC           :: ICtrDay
 PUBLIC :: Init_Class_Module_TimeManager
 PUBLIC :: SetTimeControl
 PUBLIC :: TimeIncrementSeg
 PUBLIC :: GetRec2ReadWrite
CONTAINS
 SUBROUTINE Init_Class_Module_TimeManager()
  IMPLICIT NONE
  idatei(1) = 00 
  idatei(2) = 18 
  idatei(3) = 04 
  idatei(4) = 2023 
  
  idatec(1) = 00 
  idatec(2) = 18 
  idatec(3) = 04 
  idatec(4) = 2023

  idatef(1) = 01
  idatef(2) = 24 
  idatef(3) = 05 
  idatef(4) = 2023


  DHFCT=24
  !     DT                    90
  !C * ----- < 1.0  =  100* ------ 
  !     DX                   25000  
  dt_step=240.0
  ICtrDay=0
 END SUBROUTINE Init_Class_Module_TimeManager

 FUNCTION TimeIncrementSeg(idate,idatec, ifday, tod, ihr, iday, mon, iyr,&
                            ktm,kt,ktp,ahour,bhour,maxtim,dt_step,jdt)  RESULT (ok)
    !
    !
    !==========================================================================
    ! idate(4).......date of current data
    !                idate(1)....hour(00/12)
    !                idate(2)....day of month
    !                idate(3)....month
    !                idate(4)....year
    !    ifday.......model forecast day
    !    tod.........todx=tod+swint*f3600, model forecast time of
    !                day in seconds
    !                swint....sw subr. call interval in hours
    !                swint has to be less than or equal to trint
    !                              and mod(trint,swint)=0
    !                f3600=3.6e3
    !    ihr.........hour(00/12)
    !    iday........day of month
    !    mon.........month
    !    iyr.........year
    !    yrl.........length of year in days
    !    monl(12)....length of each month in days
    !==========================================================================
    !

    INTEGER      , INTENT(in   ) :: idate(4)
    INTEGER      , INTENT(inout) :: idatec(4)
    INTEGER      , INTENT(in   ) :: ifday
    REAL(KIND=r8), INTENT(in   ) :: tod
    INTEGER      , INTENT(out  ) :: ihr
    INTEGER      , INTENT(out  ) :: iday
    INTEGER      , INTENT(out  ) :: mon
    INTEGER      , INTENT(out  ) :: iyr
    INTEGER      , INTENT(inout) :: ktm
    INTEGER      , INTENT(inout) :: kt  
    INTEGER      , INTENT(inout) :: ktp 
    REAL(KIND=r8), INTENT(inout) :: ahour
    REAL(KIND=r8), INTENT(inout) :: bhour
    INTEGER      , INTENT(in   ) :: maxtim
    REAL(KIND=r8), INTENT(in   ) :: dt_step
    INTEGER      , INTENT(in   ) :: jdt
    INTEGER          :: ioptn     
    INTEGER          :: day
    INTEGER          :: month
    INTEGER          :: year
    INTEGER          :: LenYearbyDay
    INTEGER          :: kday
    INTEGER          :: idaymn
    REAL(KIND=r8)    :: ctim
    REAL(KIND=r8)    :: hrmodl
    INTEGER          :: monl(12)
    INTEGER          :: ok

    REAL(KIND=r8), PARAMETER :: yrl =   365.2500
    REAL(KIND=r8), PARAMETER ::  ep = .015625
    DATA MONL/31,28,31,30,31,30,&
         31,31,30,31,30,31/

    ktm=kt
    ctim=tod+idate(1)*3600.0_r8
    bhour=ahour
    IF (ctim >= 86400.e0_r8) THEN
       kday=1
       ctim=ctim-86400.e0_r8
    ELSE
       kday=0
    END IF
    !
    !     adjust time to reduce round off error in divsion
    !
    iday = idate(2) + ifday + kday
    hrmodl = (ctim+ep)/3600.0_r8
    ihr = INT(hrmodl,KIND=i4)
    mon = idate(3)
    iyr = idate(4)
    DO
       ioptn=1    
       day  =31
       month  =12
       year=iyr
       CALL calndr (ioptn,  day, month, year, LenYearbyDay)
       idaymn = monl(mon)
       IF (LenYearbyDay == 366 .AND. mon == 2)idaymn=29
       IF (iday <= idaymn) exit
       iday = iday - idaymn
       mon = mon + 1
       IF (mon < 13) CYCLE
       mon = 1
       iyr = iyr + 1
    END DO
    idatec(1)=ihr
    idatec(2)=iday
    idatec(3)=mon
    idatec(4)=iyr

    ahour=(ifday*24.0e0_r8)+(tod/3.6e3_r8)
    kt   =INT(ahour-(1.0e-2_r8))
    ktp  =INT(ahour+(dt_step/3.6e3_r8)-(1.0e-2_r8))
    IF(jdt.EQ.maxtim) THEN
      ktm=kt
    END IF 

    ok=0

  END FUNCTION TimeIncrementSeg

 FUNCTION GetRec2ReadWrite(idate ,idatec)   RESULT (rec) 
  IMPLICIT NONE
  INTEGER, INTENT(IN     ) :: idate (4)
  INTEGER, INTENT(IN     ) :: idatec(4)
  INTEGER                  :: rec
  INTEGER                  :: yi
  INTEGER                  :: mi
  INTEGER                  :: di
  INTEGER                  :: hi
  INTEGER                  :: yc
  INTEGER                  :: mc
  INTEGER                  :: dc
  INTEGER                  :: hc
  REAL(KIND=r8)            :: xday
  REAL(KIND=r8)            :: xday2
  INTEGER                  :: LenYearbyDay1
  INTEGER                  :: LenYearbyDay2
  REAL(KIND=r8)            :: datehr
  REAL(KIND=r8)            :: datehf

  hi = idate (1)
  di = idate (2)
  mi = idate (3)
  yi = idate (4)
  
  hc = idatec(1)
  dc = idatec(2)
  mc = idatec(3)
  yc = idatec(4)

  rec=0
  CALL jull_options(yi,mi,di,hi,xday,LenYearbyDay1)
  datehr=(yi)*365.25e0_r8 + (xday)
  xday2=xday
  CALL jull_options(yc,mc,dc,hc,xday,LenYearbyDay2)
  datehf=(yc)*365.25e0_r8 + (xday)
  rec=nint(((datehf-datehr))*24)+1
 END FUNCTION GetRec2ReadWrite 
  
 FUNCTION SetTimeControl(idate ,idatef)   RESULT (maxtim)
  IMPLICIT NONE

  INTEGER, INTENT(IN     ) :: idate (4)
  INTEGER, INTENT(IN     ) :: idatef(4)
  INTEGER                  :: maxtim
  INTEGER                  :: yi
  INTEGER                  :: mi
  INTEGER                  :: di
  INTEGER                  :: hi
  INTEGER                  :: yf
  INTEGER                  :: mf
  INTEGER                  :: df
  INTEGER                  :: hf
  REAL(KIND=r8)            :: xday
  INTEGER                  :: LenYearbyDay
  REAL(KIND=r8)            :: datehr
  REAL(KIND=r8)            :: datehf
  INTEGER                  :: nday
  INTEGER                  :: md(12)
  INTEGER                  :: ntstep 
  REAL(KIND=r8)            :: dh
  REAL(KIND=r8)            :: nts
  REAL(KIND=r8)            :: mhf
  REAL(KIND=r8)            :: chk
  REAL(KIND=r8)            :: ntstepmax

  hi = idate (1)
  di = idate (2)
  mi = idate (3)
  yi = idate (4)
  
  hf = idatef(1)
  df = idatef(2)
  mf = idatef(3)
  yf = idatef(4)

  CALL jull_options(yi,mi,di,hi,xday,LenYearbyDay)
  datehr=yi+(xday/365.25e0_r8)
  CALL jull_options(yf,mf,df,hf,xday,LenYearbyDay)
  datehf=yf+(xday/365.25e0_r8)

  nday=0
    IF(yi == yf .AND. mi==mf .AND. di==df) THEN
       nday=0
    ELSE
       CALL jull_options(yi,mi,di,hi,xday,LenYearbyDay)
       DO WHILE (datehr <= datehf)
          nday=nday+1
          IF ( LenYearbyDay == 366 )THEN
             md =(/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
          ELSE
             md =(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
          END IF
          di=di+1
          IF( di > md(mi) )THEN
             di=1
             mi=mi+1
             IF ( mi > 12 ) THEN
                mi=1
                yi=yi+1
             END IF
          END IF
          CALL jull_options(yi,mi,di,hi,xday,LenYearbyDay)
          datehr=yi+(xday/365.25e0_r8)
       END DO
    END IF
    ntstep   =INT(REAL((nday)*86400)/dt_step)

   IF ( dhfct /= 0 ) THEN
      IF ( hi /= hf ) THEN
         dh =hf-hi
         nts=dh*3600/dt_step
         mhf=dh/dhfct
         chk=mhf*dhfct
         IF ( chk /= dh ) THEN
            WRITE(nfprt,*) 'Wrong Request for the Hour in datef =', yf,mf,df,hf
            WRITE(nfprt,*) 'Difference of Hours in datei = ',yi,mi,di,hi, 'and '
            WRITE(nfprt,*) 'datef is Not Compatible With dhfct =' ,dhfct
            STOP
         END IF
         ntstep=ntstep+INT(nts)
      END IF
   END IF
   maxtim=ntstep
   ntstepmax=INT(REAL(51*366*86400)/dt_step)
   cth0 =dhfct

    IF( ntstep > ntstepmax ) THEN
       WRITE(nfprt,*) 'nstep = ',ntstep,' is greater than ntstepmax = ',ntstepmax
       STOP
    END IF
    dct=cth0

  END FUNCTION SetTimeControl


  SUBROUTINE jull_options(yi,mi,di,hi,xday,LenYearbyDay)

    INTEGER, INTENT(IN   ) :: yi
    INTEGER, INTENT(IN   ) :: mi
    INTEGER, INTENT(IN   ) :: di
    INTEGER, INTENT(IN   ) :: hi
    REAL(KIND=r8)   , INTENT(OUT  ) :: xday
    INTEGER         , INTENT(OUT  ) :: LenYearbyDay
    REAL(KIND=r8)                   :: tod
    REAL(KIND=r8)                   :: yrl
    !INTEGER                :: monl(12)
    !INTEGER                :: monday(12)
    INTEGER                :: idayct2,day, mon, year,ioptn
    REAL(KIND=r8)   , PARAMETER     :: f3600=3.6e3_r8
    ioptn=1    
    day=di
    mon=mi
    year=yi
    CALL calndr (ioptn,  day, mon, year, idayct2)
    tod=0.0_r8
    yrl=365.25e0_r8
    !MONL    =   (/31,28,31,30,31,30,31,31,30,31,30,31/)
    !
    !     id is now assumed to be the current date and hour
    !
    !monday(1)=0
    !DO m=2,12
    !   monday(m)=monday(m-1)+monl(m-1)
    !END DO
    xday=hi*f3600
    xday=xday+MOD(tod,f3600)
    IF(hi > 0 ) THEN
      xday=(idayct2-1) + xday/86400.0_r8
    ELSE
      xday=idayct2 + xday/86400.0_r8
    END IF
    ioptn=1    
    day=31
    mon=12
    year=yi+3
    CALL calndr (ioptn,  day, mon, year, LenYearbyDay)
    IF(LenYearbyDay /= 366)THEN
       xday=xday-MOD(yi+3,4)*0.25_r8
    ELSE
       xday=xday-(0*0.25_r8)
    END IF
    ioptn=1    
    day=31
    mon=12
    year=yi
    CALL calndr (ioptn,  day, mon, year, LenYearbyDay)
    !IF(LenYearbyDay == 366 .AND.mi.GT.2)xday=xday+1.0e0_r8
    xday= MOD(xday-1.0_r8,yrl)

  END SUBROUTINE jull_options


  SUBROUTINE calndr(ioptn,iday,month,iyear,idayct)

    ! CALNDR = CALeNDaR conversions, version 1.0

    IMPLICIT NONE

    ! specify the desired calendar conversion option.
    ! in order to return the julian day number, compatible with function idaywk from above,
    ! we choose option 3
    ! (tested with dates: Feb, 23 2010 -> idaywk = Tue
    !                               Dec, 24 2009 -> idaywk = Thu
    !                               Oct, 15 1582  -> idaywk = Fri ...which all look o.k. )
    INTEGER, INTENT(in) :: ioptn 

    ! Input/Output variables
    INTEGER, INTENT(inout) :: iday,month,iyear,idayct

    !----------
    !
    ! The subroutine calndr() performs calendar calculations using either
    ! the standard Gregorian calendar or the old Julian calendar.
    ! This subroutine extends the definitions of these calendar systems
    ! to any arbitrary year.  The algorithms in this subroutine
    ! will work with any date in the past or future,
    ! but overflows will occur if the numbers are sufficiently large.
    ! For a computer using a 32-bit integer, this routine can handle
    ! any date between roughly 5.8 million BC and 5.8 million AD
    ! without experiencing overflow during calculations.
    !
    ! No external functions or subroutines are called.
    !
    !----------
    !
    ! Input/output arguments for subroutine CALNDR()
    !
    ! "ioptn" is the desired calendar conversion option explained below.
    ! Positive option values use the standard modern Gregorian calendar.
    ! Negative option values use the old Julian calendar which was the
    ! standard in Europe from its institution by Julius Caesar in 45 BC
    ! until at least 4 October 1582.  The Gregorian and Julian calendars
    ! are explained further below.
    !
    ! (iday,month,iyear) is a calendar date where "iday" is the day of
    ! the month, "month" is 1 for January, 2 for February, etc.,
    ! and "iyear" is the year.  If the year is 1968 AD, enter iyear=1968,
    ! since iyear=68 would refer to 68 AD.
    ! For BC years, iyear should be negative, so 45 BC would be iyear=-45.
    ! By convention, there is no year 0 under the BC/AD year numbering
    ! scheme.  That is, years proceed as 2 BC, 1 BC, 1 AD, 2 AD, etc.,
    ! without including 0. The subroutine calndr() will print an error message
    ! and stop if you specify iyear = 0.
    !
    ! "idayct" is a day count.  It is either the day number during the
    ! specified year or the Julian Day number, depending on the value
    ! of ioptn.  By day number during the specified year, we mean
    ! idayct=1 on 1 January, idayct=32 on 1 February, etc., to idayct=365
    ! or 366 on 31 December, depending on whether the specified year
    ! is a leap year.
    !
    ! The values of input variables are not changed by this subroutine.
    !
    !
    ! ALLOWABLE VALUES FOR "IOPTN" and the conversions they invoke.
    ! Positive option values ( 1 to  5) use the standard Gregorian calendar.
    ! Negative option values (-1 to -5) use the old      Julian    calendar.
    !
    ! Absolute
    !  value
    ! of ioptn   Input variable(s)     Output variable(s)
    !
    !    1       iday,month,iyear      idayct
    ! Given a calendar date (iday,month,iyear), compute the day number
    ! (idayct) during the year, where 1 January is day number 1 and
    ! 31 December is day number 365 or 366, depending on whether it is
    ! a leap year.
    !
    !    2       idayct,iyear          iday,month
    ! Given the day number of the year (idayct) and the year (iyear),
    ! compute the day of the month (iday) and the month (month).
    !
    !    3       iday,month,iyear      idayct
    ! Given a calendar date (iday,month,iyear), compute the Julian Day
    ! number (idayct) that starts at noon of the calendar date specified.
    !
    !    4       idayct                iday,month,iyear
    ! Given the Julian Day number (idayct) that starts at noon,
    ! compute the corresponding calendar date (iday,month,iyear).
    !
    !    5       idayct                iday,month,iyear
    ! Given the Julian Day number (idayct) that starts at noon,
    ! compute the corresponding day number for the year (iday)
    ! and year (iyear).  On return from calndr(), "month" will always
    ! be set equal to 1 when ioptn=5.
    !
    ! No inverse function is needed for ioptn=5 because it is
    ! available through option 3.  One simply calls calndr() with:
    ! ioptn = 3,
    ! iday  = day number of the year instead of day of the month,
    ! month = 1, and
    ! iyear = whatever the desired year is.
    !
    !----------
    !
    ! EXAMPLES
    ! The first 6 examples are for the standard Gregorian calendar.
    ! All the examples deal with 15 October 1582, which was the first day
    ! of the Gregorian calendar.  15 October is the 288-th day of the year.
    ! Julian Day number 2299161 began at noon on 15 October 1582.
    !
    ! Find the day number during the year on 15 October 1582
    !     ioptn = 1
    !     call calndr (ioptn, 15, 10, 1582,  idayct)
    ! calndr() should return idayct=288
    !
    ! Find the day of the month and month for day 288 in year 1582.
    !     ioptn = 2
    !     call calndr (ioptn, iday, month, 1582, 288)
    ! calndr() should return iday=15 and month=10.
    !
    ! Find the Julian Day number for 15 October 1582.
    !     ioptn = 3
    !     call calndr (ioptn, 15, 10, 1582, julian)
    ! calndr() should return julian=2299161
    !
    ! Find the Julian Day number for day 288 during 1582 AD.
    ! When the input is day number of the year, one should specify month=1
    !     ioptn = 3
    !     call calndr (ioptn, 288, 1, 1582, julian)
    ! calndr() should return dayct=2299161
    !
    ! Find the date for Julian Day number 2299161.
    !     ioptn = 4
    !     call calndr (ioptn, iday, month, iyear, 2299161)
    ! calndr() should return iday=15, month=10, and iyear=1582
    !
    ! Find the day number during the year (iday) and year
    ! for Julian Day number 2299161.
    !     ioptn = 5
    !     call calndr (ioptn, iday, month, iyear, 2299161)
    ! calndr() should return iday=288, month = 1, iyear=1582
    !
    ! Given 15 October 1582 under the Gregorian calendar,
    ! find the date (idayJ,imonthJ,iyearJ) under the Julian calendar.
    ! To do this, we call calndr() twice, using the Julian Day number
    ! as the intermediate value.
    !     call calndr ( 3, 15,        10, 1582,    julian)
    !     call calndr (-4, idayJ, monthJ, iyearJ,  julian)
    ! The first call to calndr() should return julian=2299161, and
    ! the second should return idayJ=5, monthJ=10, iyearJ=1582
    !
    !----------
    !
    ! BASIC CALENDAR INFORMATION
    !
    ! The Julian calendar was instituted by Julius Caesar in 45 BC.
    ! Every fourth year is a leap year in which February has 29 days.
    ! That is, the Julian calendar assumes that the year is exactly
    ! 365.25 days long.  Actually, the year is not quite this long.
    ! The modern Gregorian calendar remedies this by omitting leap years
    ! in years divisible by 100 except when the year is divisible by 400.
    ! Thus, 1700, 1800, and 1900 are leap years under the Julian calendar
    ! but not under the Gregorian calendar.  The years 1600 and 2000 are
    ! leap years under both the Julian and the Gregorian calendars.
    ! Other years divisible by 4 are leap years under both calendars,
    ! such as 1992, 1996, 2004, 2008, 2012, etc.  For BC years, we recall
    ! that year 0 was omitted, so 1 BC, 5 BC, 9 BC, 13 BC, etc., and 401 BC,
    ! 801 BC, 1201 BC, etc., are leap years under both calendars, while
    ! 101 BC, 201 BC, 301 BC, 501 BC, 601 BC, 701 BC, 901 BC, 1001 BC,
    ! 1101 BC, etc., are leap years under the Julian calendar but not
    ! the Gregorian calendar.
    !
    ! The Gregorian calendar is named after Pope Gregory XIII.  He declared
    ! that the last day of the old Julian calendar would be Thursday,
    ! 4 October 1582 and that the following day, Friday, would be reckoned
    ! under the new calendar as 15 October 1582.  The jump of 10 days was
    ! included to make 21 March closer to the spring equinox.
    !
    ! Only a few Catholic countries (Italy, Poland, Portugal, and Spain)
    ! switched to the Gregorian calendar on the day after 4 October 1582.
    ! It took other countries months to centuries to change to the
    ! Gregorian calendar.  For example, England's first day under the
    ! Gregorian calendar was 14 September 1752.  The same date applied to
    ! the entire British empire, including America.  Japan, Russia, and many
    ! eastern European countries did not change to the Gregorian calendar
    ! until the 20th century.  The last country to change was Turkey,
    ! which began using the Gregorian calendar on 1 January 1927.
    !
    ! Therefore, between the years 1582 and 1926 AD, you must know
    ! the country in which an event was dated to interpret the date
    ! correctly.  In Sweden, there was even a year (1712) when February
    ! had 30 days.  Consult a book on calendars for more details
    ! about when various countries changed their calendars.
    !
    ! DAY NUMBER DURING THE YEAR
    ! The day number during the year is simply a counter equal to 1 on
    ! 1 January, 32 on 1 February, etc., through 365 or 366 on 31 December,
    ! depending on whether the year is a leap year.  Sometimes this is
    ! called the Julian Day, but that term is better reserved for the
    ! day counter explained below.
    !
    ! JULIAN DAY NUMBER
    ! The Julian Day numbering system was designed by Joseph Scaliger
    ! in 1582 to remove ambiguity caused by varying calendar systems.
    ! The name "Julian Day" was chosen to honor Scaliger's father,
    ! Julius Caesar Scaliger (1484-1558), an Italian scholar and physician
    ! who lived in France.  Because Julian Day numbering was especially
    ! designed for astronomers, Julian Days begin at noon so that the day
    ! counter does not change in the middle of an astronomer's observing
    ! period.  Julian Day 0 began at noon on 1 January 4713 BC under the
    ! Julian calendar.  A modern reference point is that 23 May 1968
    ! (Gregorian calendar) was Julian Day 2,440,000.
    !
    ! JULIAN DAY NUMBER EXAMPLES
    !
    ! The table below shows a few Julian Day numbers and their corresponding
    ! dates, depending on which calendar is used.  A negative 'iyear' refers
    ! to BC (Before Christ).
    !
    !                     Julian Day under calendar:
    ! iday  month   iyear     Gregorian   Julian
    !  24     11   -4714            0        -38
    !   1      1   -4713           38          0
    !   1      1       1      1721426    1721424
    !   4     10    1582      2299150    2299160
    !  15     10    1582      2299161    2299171
    !   1      3    1600      2305508    2305518
    !  23      5    1968      2440000    2440013
    !   5      7    1998      2451000    2451013
    !   1      3    2000      2451605    2451618
    !   1      1    2001      2451911    2451924
    !
    ! From this table, we can see that the 10 day difference between the
    ! two calendars in 1582 grew to 13 days by 1 March 1900, since 1900 was
    ! a leap year under the Julian calendar but not under the Gregorian
    ! calendar.  The gap will widen to 14 days after 1 March 2100 for the
    ! same reason.
    !
    !----------
    !
    ! PORTABILITY
    !
    ! This subroutine is written in standard Fortran 90.
    ! It calls no external functions or subroutines and should run
    ! without problem on any computer having a 32-bit word or longer.
    !
    !----------
    !
    ! ALGORITHM
    !
    ! The goal in coding calndr() was clear, clean code, not efficiency.
    ! Calendar calculations usually take a trivial fraction of the time
    ! in any program in which dates conversions are involved.
    ! Data analysis usually takes the most time.
    !
    ! Standard algorithms are followed in this subroutine.  Internal to
    ! this subroutine, we use a year counter "jyear" such that
    !  jyear=iyear   when iyear is positive
    !       =iyear+1 when iyear is negative.
    ! Thus, jyear does not experience a 1 year jump like iyear does
    ! when going from BC to AD.  Specifically, jyear = 0 when iyear=-1,
    ! i.e., when the year is 1 BC.
    !
    ! For simplicity in dealing with February, inside this subroutine,
    ! we let the year begin on 1 March so that the adjustable month,
    ! February is the last month of the year.
    ! It is clear that the calendar used to work this way because the
    ! months September, October, November, and December refer to
    ! 7, 8, 9, and 10.  For consistency, jyear is incremented on 1 March
    ! rather than on 1 January.  Of course, everything is adjusted back to
    ! standard practice of years beginning on 1 January before answers
    ! are returned to the routine that calls calndr().
    !
    ! Lastly, we use a trick to calculate the number of days from 1 March
    ! until the end of the month that precedes the specified month.
    ! That number of days is int(30.6001*(month+1))-122,
    ! where 30.6001 is used to avoid the possibility of round-off and
    ! truncation error.  For example, if 30.6 were used instead,
    ! 30.6*5 should be 153, but round-off error could make it 152.99999,
    ! which would then truncated to 152, causing an error of 1 day.
    !
    ! Algorithm reference:
    ! Dershowitz, Nachum and Edward M. Reingold, 1990: Calendrical
    ! Calculations.  Software-Practice and Experience, vol. 20, number 9
    ! (September 1990), pp. 899-928.
    !
    ! Copyright (C) 1999 Jon Ahlquist.
    ! Issued under the second GNU General Public License.
    ! See www.gnu.org for details.
    ! This program is distributed in the hope that it will be useful,
    ! but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    ! If you find any errors, please notify:
    ! Jon Ahlquist
    ! Dept of Meteorology
    ! Florida State University
    ! Tallahassee, FL 32306-4520
    ! 15 March 1999.
    !
    !-----

    ! converted to Fortran90 by Dimitri Komatitsch,
    ! University of Pau, France, January 2008.

    ! Declare internal variables.
    INTEGER jdref, jmonth, jyear, leap, n1yr, n4yr, n100yr, n400yr, ndays, ndy400, ndy100, nyrs, yr400, yrref
    !
    ! Explanation of all internal variables.
    ! jdref   Julian Day on which 1 March begins in the reference year.
    ! jmonth  Month counter which equals month+1 if month > 2
    !          or month+13 if month <= 2.
    ! jyear   Year index,  jyear=iyear if iyear > 0, jyear=iyear+1
    !            if iyear < 0.  Thus, jyear does not skip year 0
    !            like iyear does between BC and AD years.
    ! leap    =1 if the year is a leap year,  = 0 if not.
    ! n1yr    Number of complete individual years between iyear and
    !            the reference year after all 4, 100,
    !            and 400 year periods have been removed.
    ! n4yr    Number of complete 4 year cycles between iyear and
    !            the reference year after all 100 and 400 year periods
    !            have been removed.
    ! n100yr  Number of complete 100 year periods between iyear and
    !            the reference year after all 400 year periods
    !            have been removed.
    ! n400yr  Number of complete 400 year periods between iyear and
    !            the reference year.
    ! ndays   Number of days since 1 March during iyear.  (In intermediate
    !            steps, it holds other day counts as well.)
    ! ndy400  Number of days in 400 years.  Under the Gregorian calendar,
    !            this is 400*365 + 100 - 3 = 146097.  Under the Julian
    !            calendar, this is 400*365 + 100 = 146100.
    ! ndy100  Number of days in 100 years,  Under the Gregorian calendar,
    !            this is 100*365 + 24 = 36524.   Under the Julian calendar,
    !            this is 100*365 + 25 = 36525.
    ! nyrs    Number of years from the beginning of yr400
    !              to the beginning of jyear.  (Used for option +/-3).
    ! yr400   The largest multiple of 400 years that is <= jyear.
    !
    !
    !----------------------------------------------------------------
    ! Do preparation work.
    !
    ! Look for out-of-range option values.
    IF ((ioptn == 0) .OR. (ABS(ioptn) >= 6)) THEN
       WRITE(*,*)'For calndr(), you specified ioptn = ', ioptn
       WRITE(*,*) 'Allowable values are 1 to 5 for the Gregorian calendar'
       WRITE(*,*) 'and -1 to -5 for the Julian calendar.'
       STOP
    ENDIF
    !
    ! Options 1-3 have "iyear" as an input value.
    ! Internally, we use variable "jyear" that does not have a jump
    ! from -1 (for 1 BC) to +1 (for 1 AD).
    IF (ABS(ioptn) <= 3) THEN
       IF (iyear > 0) THEN
          jyear = iyear
       ELSE IF (iyear == 0) THEN
          WRITE(*,*) 'For calndr(), you specified the nonexistent year 0'
          STOP
       ELSE
          jyear = iyear + 1
       ENDIF
       !
       !        Set "leap" equal to 0 if "jyear" is not a leap year
       !        and equal to 1 if it is a leap year.
       leap = 0
       IF ((jyear/4)*4 == jyear) THEN
          leap = 1
       ENDIF
       IF ((ioptn > 0) .AND. ((jyear/100)*100 == jyear) .AND. ((jyear/400)*400 /= jyear)) THEN
          leap = 0
       ENDIF
    ENDIF
    !
    ! Options 3-5 involve Julian Day numbers, which need a reference year
    ! and the Julian Days that began at noon on 1 March of the reference
    ! year under the Gregorian and Julian calendars.  Any year for which
    ! "jyear" is divisible by 400 can be used as a reference year.
    ! We chose 1600 AD as the reference year because it is the closest
    ! multiple of 400 to the institution of the Gregorian calendar, making
    ! it relatively easy to compute the Julian Day for 1 March 1600
    ! given that, on 15 October 1582 under the Gregorian calendar,
    ! the Julian Day was 2299161.  Similarly, we need to do the same
    ! calculation for the Julian calendar.  We can compute this Julian
    ! Day knowing that on 4 October 1582 under the Julian calendar,
    ! the Julian Day number was 2299160.  The details of these calculations
    ! is next.
    !    From 15 October until 1 March, the number of days is the remainder
    ! of October plus the days in November, December, January, and February:
    ! 17+30+31+31+28 = 137, so 1 March 1583 under the Gregorian calendar
    ! was Julian Day 2,299,298.  Because of the 10 day jump ahead at the
    ! switch from the Julian calendar to the Gregorian calendar, 1 March
    ! 1583 under the Julian calendar was Julian Day 2,299,308.  Making use
    ! of the rules for the two calendar systems, 1 March 1600 was Julian
    ! Day 2,299,298 + (1600-1583)*365 + 5 (due to leap years) =
    ! 2,305,508 under the Gregorian calendar and day 2,305,518 under the
    ! Julian calendar.
    !    We also set the number of days in 400 years and 100 years.
    ! For reference, 400 years is 146097 days under the Gregorian calendar
    ! and 146100 days under the Julian calendar.  100 years is 36524 days
    ! under the Gregorian calendar and 36525 days under the Julian calendar.
    IF (ABS(ioptn) >= 3) THEN
       !
       !        Julian calendar values.
       yrref  =    1600
       jdref  = 2305518
       !               = Julian Day reference value for the day that begins
       !                 at noon on 1 March of the reference year "yrref".
       ndy400 = 400*365 + 100
       ndy100 = 100*365 +  25
       !
       !        Adjust for Gregorian calendar values.
       IF (ioptn > 0) THEN
          jdref  = jdref  - 10
          ndy400 = ndy400 -  3
          ndy100 = ndy100 -  1
       ENDIF
    ENDIF
    !
    !----------------------------------------------------------------
    ! OPTIONS -1 and +1:
    ! Given a calendar date (iday,month,iyear), compute the day number
    ! of the year (idayct), where 1 January is day number 1 and 31 December
    ! is day number 365 or 366, depending on whether it is a leap year.
    IF (ABS(ioptn) == 1) THEN
       !
       !     Compute the day number during the year.
       IF (month <= 2) THEN
          idayct = iday + (month-1)*31
       ELSE
          idayct = iday + INT(30.6001 * (month+1)) - 63 + leap
       ENDIF
       !
       !----------------------------------------------------------------
       ! OPTIONS -2 and +2:
       ! Given the day number of the year (idayct) and the year (iyear),
       ! compute the day of the month (iday) and the month (month).
    ELSE IF (ABS(ioptn) == 2) THEN
       !
       IF (idayct < 60+leap) THEN
          month  = (idayct-1)/31
          iday   = idayct - month*31
          month  = month + 1
       ELSE
          ndays  = idayct - (60+leap)
          !               = number of days past 1 March of the current year.
          jmonth = (10*(ndays+31))/306 + 3
          !               = month counter, =4 for March, =5 for April, etc.
          iday   = (ndays+123) - INT(30.6001*jmonth)
          month  = jmonth - 1
       ENDIF
       !
       !----------------------------------------------------------------
       ! OPTIONS -3 and +3:
       ! Given a calendar date (iday,month,iyear), compute the Julian Day
       ! number (idayct) that starts at noon.
    ELSE IF (ABS(ioptn) == 3) THEN
       !
       !     Shift to a system where the year starts on 1 March, so January
       !     and February belong to the preceding year.
       !     Define jmonth=4 for March, =5 for April, ..., =15 for February.
       IF (month <= 2) THEN
          jyear  = jyear -  1
          jmonth = month + 13
       ELSE
          jmonth = month +  1
       ENDIF
       !
       !     Find the closest multiple of 400 years that is <= jyear.
       yr400 = (jyear/400)*400
       !           = multiple of 400 years at or less than jyear.
       IF (jyear < yr400) THEN
          yr400 = yr400 - 400
       ENDIF
       !
       n400yr = (yr400 - yrref)/400
       !            = number of 400-year periods from yrref to yr400.
       nyrs   = jyear - yr400
       !            = number of years from the beginning of yr400
       !              to the beginning of jyear.
       !
       !     Compute the Julian Day number.
       idayct = iday + INT(30.6001*jmonth) - 123 + 365*nyrs + nyrs/4 &
            + jdref + n400yr*ndy400
       !
       !     If we are using the Gregorian calendar, we must not count
       !     every 100-th year as a leap year.  nyrs is less than 400 years,
       !     so we do not need to consider the leap year that would occur if
       !     nyrs were divisible by 400, i.e., we do not add nyrs/400.
       IF (ioptn > 0) THEN
          idayct = idayct - nyrs/100
       ENDIF
       !
       !----------------------------------------------------------------
       ! OPTIONS -5, -4, +4, and +5:
       ! Given the Julian Day number (idayct) that starts at noon,
       ! compute the corresponding calendar date (iday,month,iyear)
       ! (abs(ioptn)=4) or day number during the year (abs(ioptn)=5).
    ELSE
       !
       !     Create a new reference date which begins on the nearest
       !     400-year cycle less than or equal to the Julian Day for 1 March
       !     in the year in which the given Julian Day number (idayct) occurs.
       ndays  = idayct - jdref
       n400yr = ndays / ndy400
       !            = integral number of 400-year periods separating
       !              idayct and the reference date, jdref.
       jdref  = jdref + n400yr*ndy400
       IF (jdref > idayct) THEN
          n400yr = n400yr - 1
          jdref  = jdref  - ndy400
       ENDIF
       !
       ndays  = idayct - jdref
       !            = number from the reference date to idayct.
       !
       n100yr = MIN(ndays/ndy100, 3)
       !            = number of complete 100-year periods
       !              from the reference year to the current year.
       !              The min() function is necessary to avoid n100yr=4
       !              on 29 February of the last year in the 400-year cycle.
       !
       ndays  = ndays - n100yr*ndy100
       !            = remainder after removing an integral number of
       !              100-year periods.
       !
       n4yr   = ndays / 1461
       !            = number of complete 4-year periods in the current century.
       !              4 years consists of 4*365 + 1 = 1461 days.
       !
       ndays  = ndays - n4yr*1461
       !            = remainder after removing an integral number
       !              of 4-year periods.
       !
       n1yr   = MIN(ndays/365, 3)
       !            = number of complete years since the last leap year.
       !              The min() function is necessary to avoid n1yr=4
       !              when the date is 29 February on a leap year,
       !              in which case ndays=1460, and 1460/365 = 4.
       !
       ndays  = ndays - 365*n1yr
       !            = number of days so far in the current year,
       !              where ndays = 0 on 1 March.
       !
       iyear  = n1yr + 4*n4yr + 100*n100yr + 400*n400yr + yrref
       !            = year, as counted in the standard way,
       !              but relative to 1 March.
       !
       ! At this point, we need to separate ioptn=abs(4), which seeks a
       ! calendar date, and ioptn=abs(5), which seeks the day number during
       ! the year.  First compute the calendar date if desired (abs(ioptn)=4).
       IF (ABS(ioptn) == 4) THEN
          jmonth = (10*(ndays+31))/306 + 3
          !               = offset month counter.  jmonth=4 for March, =13 for
          !                 December, =14 for January, =15 for February.
          iday   = (ndays+123) - INT(30.6001*jmonth)
          !               = day of the month, starting with 1 on the first day
          !                 of the month.
          !
          !        Now adjust for the fact that the year actually begins
          !        on 1 January.
          IF (jmonth <= 13) THEN
             month = jmonth - 1
          ELSE
             month = jmonth - 13
             iyear = iyear + 1
          ENDIF
          !
          ! This code handles abs(ioptn)=5, finding the day number during the year.
       ELSE
          !        ioptn=5 always returns month = 1, which we set now.
          month = 1
          !
          !        We need to determine whether this is a leap year.
          leap = 0
          IF ((jyear/4)*4 == jyear) THEN
             leap = 1
          ENDIF
          IF ((ioptn > 0) .AND. ((jyear/100)*100 == jyear) .AND. ((jyear/400)*400 /= jyear)) THEN
             leap = 0
          ENDIF
          !
          !        Now find the day number "iday".
          !        ndays is the number of days since the most recent 1 March,
          !        so ndays = 0 on 1 March.
          IF (ndays <= 305) THEN
             iday  = ndays + 60 + leap
          ELSE
             iday  = ndays - 305
             iyear = iyear + 1
          ENDIF
       ENDIF
       !
       !     Adjust the year if it is <= 0, and hence BC (Before Christ). 
       IF (iyear <= 0) THEN
          iyear = iyear - 1
       ENDIF
       !
       ! End the code for the last option, ioptn.
    ENDIF

  END SUBROUTINE calndr

END MODULE Class_Module_TimeManager
