module FPpdecol
contains
!///////////////////////////////////////////////////////////////////////
  SUBROUTINE PDECOL (T0, TOUT, DT, XBKPT, EPS, NINT, KORD, NCC,     &
       NPDE, MF, INDEX, WORK, IWORK)                                     
    use FPglobal, ONLY: NOGAUS, MAXDER
    use FPwork
    use FPuser
    save ! all variables                                              
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
! THIS IS THE JULY 12, 1977 VERSION OF PDECOL.                          
!                                                                       
! THIS PACKAGE WAS CONSTRUCTED SO AS TO CONFORM TO AS MANY ANSI-FORTRAN 
! RULES AS WAS CONVIENTLY POSSIBLE.  THE FORTRAN USED VIOLATES ANSI     
! STANDARDS IN THE TWO WAYS LISTED BELOW....                            
!                                                                       
!   1. SUBSCRIPTS OF THE GENERAL FORM C*V1 + V2 + V3 ARE USED           
!      (POSSIBLY IN A PERMUTED ORDER), WHERE C IS AN INTEGER CONSTANT   
!      AND V1, V2, AND V3 ARE INTEGER VARIABLES.                        
!                                                                       
!   2. ARRAY NAMES APPEAR SINGLY IN DATA STATEMENTS IN THE ROUTINES     
!      BSPLVN AND CSTCOL.                                               
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
!-----------------------------------------------------------------------
! PDECOL IS THE DRIVER ROUTINE FOR A SOPHISTICATED PACKAGE OF           
! SUBROUTINES WHICH IS DESIGNED TO SOLVE THE GENERAL SYSTEM OF          
! NPDE NONLINEAR PARTIAL DIFFERENTIAL EQUATIONS OF AT MOST SECOND       
! ORDER ON THE INTERVAL (XLEFT,XRIGHT) FOR T .GT. T0 WHICH IS OF THE    
! FORM....                                                              
!                                                                       
!      DU/DT  =  F( T, X, U, UX, UXX )                                  
!                                                                       
! WHERE                                                                 
!                                                                       
!          U  =  (  U(1),  U(2), ... ,  U(NPDE) )                       
!         UX  =  ( UX(1), UX(2), ... , UX(NPDE) )                       
!        UXX  =  (UXX(1),UXX(2), ... ,UXX(NPDE) ) .                     
!                                                                       
! EACH U(K) IS A FUNCTION OF THE SCALAR QUANTITIES T AND X.             
! UX(K) REPRESENTS THE FIRST PARTIAL DERIVATIVE OF U(K) WITH RESPECT    
! TO THE VARIABLE X,  UXX(K) REPRESENTS THE SECOND PARTIAL DERIVATIVE   
! OF U(K) WITH RESPECT TO THE VARIABLE X, AND DU/DT IS THE VECTOR OF    
! PARTIAL DERIVATIVES OF U WITH RESPECT TO THE TIME VARIABLE T.         
! F  REPRESENTS AN ARBITRARY VECTOR VALUED FUNCTION WHOSE NPDE          
! COMPONENTS DEFINE THE RESPECTIVE PARTIAL DIFFERENTIAL EQUATIONS OF    
! THE PDE SYSTEM.  SEE SUBROUTINE F DESCRIPTION BELOW.                  
!                                                                       
! BOUNDARY CONDITIONS                                                   
!                                                                       
!   DEPENDING ON THE TYPE OF PDE(S), 0, 1, OR 2 BOUNDARY CONDITIONS     
!   ARE REQUIRED FOR EACH PDE IN THE SYSTEM. THESE ARE IMPOSED AT XLEFT 
!   AND/OR XRIGHT AND EACH MUST BE OF THE FORM....                      
!                                                                       
!        B(U,UX)  =  Z(T)                                               
!                                                                       
!   WHERE  B  AND  Z  ARE ARBITRARY VECTOR VALUED FUNCTIONS WITH        
!   NPDE COMPONENTS AND  U, UX, AND T  ARE AS ABOVE.  THESE BOUNDARY    
!   CONDITIONS MUST BE CONSISTENT WITH THE INITIAL CONDITIONS WHICH ARE 
!   DESCRIBED NEXT.                                                     
!                                                                       
! INITIAL CONDITIONS                                                    
!                                                                       
!   EACH SOLUTION COMPONENT  U(K)  IS ASSUMED TO BE A KNOWN (USER       
!   PROVIDED) FUNCTION OF  X  AT THE INITIAL TIME T = T0.  THE          
!   INITIAL CONDITION FUNCTIONS MUST BE CONSISTENT WITH THE BOUNDARY    
!   CONDITIONS ABOVE, I.E. THE INITIAL CONDITION FUNCTIONS MUST         
!   SATISFY THE BOUNDARY CONDITIONS FOR T = T0.  SEE SUBROUTINE UINIT   
!   DESCRIPTION BELOW.                                                  
!-----------------------------------------------------------------------
!                                                                       
! REQUIRED USER SUPPLIED SUBROUTINES                                    
!                                                                       
! THE USER IS REQUIRED TO CONSTRUCT THREE SUBPROGRAMS AND A MAIN        
! PROGRAM WHICH DEFINE THE PDE PROBLEM WHOSE SOLUTION IS TO BE          
! ATTEMPTED.  THE THREE SUBPROGRAMS ARE...                              
!                                                                       
! 1)  SUBROUTINE F( T, X, U, UX, UXX, FVAL, NPDE )                      
!     DIMENSION U(NPDE), UX(NPDE), UXX(NPDE), FVAL(NPDE)                
!        THIS ROUTINE DEFINES THE DESIRED PARTIAL DIFFERENTIAL          
!        EQUATIONS TO BE SOLVED.  THE PACKAGE PROVIDES VALUES OF THE    
!        INPUT SCALARS T AND X AND INPUT ARRAYS (LENGTH NPDE) U, UX,    
!        AND UXX, AND THE USER MUST CONSTRUCT THIS ROUTINE TO COMPUTE   
!        THE OUTPUT ARRAY FVAL (LENGTH NPDE) WHICH CONTAINS THE         
!        CORRESPONDING VALUES OF THE RIGHT HAND SIDES OF THE DESIRED    
!        PARTIAL DIFFERENTIAL EQUATIONS, I.E.                           
!                                                                       
!        FVAL(K) = THE VALUE OF THE RIGHT HAND SIDE OF THE K-TH PDE IN  
!                  THE PDE SYSTEM ABOVE, FOR K = 1 TO NPDE.             
!                                                                       
!        THE INCOMING VALUE OF THE SCALAR QUANTITY X WILL BE A          
!        COLLOCATION POINT VALUE (SEE INITAL AND COLPNT) AND THE        
!        INCOMING VALUES IN THE ARRAYS U, UX AND UXX CORRESPOND TO THIS 
!        POINT X AND TIME T.                                            
!     RETURN                                                            
!     END                                                               
!                                                                       
! 2)  SUBROUTINE BNDRY( T, X, U, UX, DBDU, DBDUX, DZDT, NPDE )          
!     DIMENSION U(NPDE), UX(NPDE), DZDT(NPDE)                           
!     DIMENSION DBDU(NPDE,NPDE), DBDUX(NPDE,NPDE)                       
!        THIS ROUTINE IS USED TO PROVIDE THE PDE PACKAGE WITH NEEDED    
!        INFORMATION ABOUT THE BOUNDARY CONDITION FUNCTIONS B AND Z     
!        ABOVE.  THE PACKAGE PROVIDES VALUES OF THE INPUT VARIABLES     
!        T, X, U, AND UX, AND THE USER IS TO DEFINE THE CORRESPONDING   
!        OUTPUT VALUES OF THE DERIVATIVES OF THE FUNCTIONS B AND Z      
!        WHERE....                                                      
!           DBDU(K,J) = PARTIAL DERIVATIVE OF THE K-TH COMPONENT OF THE 
!                       VECTOR FUNCTION B(U,UX) ABOVE WITH RESPECT TO   
!                       THE J-TH VARIABLE U(J).                         
!          DBDUX(K,J) = PARTIAL DERIVATIVE OF THE K-TH COMPONENT OF THE 
!                       VECTOR FUNCTION B(U,UX) ABOVE WITH RESPECT TO   
!                       THE J-TH VARIABLE UX(J).                        
!             DZDT(K) = DERIVATIVE OF THE K-TH COMPONENT OF THE VECTOR  
!                       FUNCTION Z(T) ABOVE WITH RESPECT TO THE         
!                       VARIABLE T.                                     
!        NOTE... THE INCOMING VALUE OF X WILL BE EITHER XLEFT OR XRIGHT.
!        IF NO BOUNDARY CONDITION IS DESIRED FOR SAY THE K-TH PDE AT    
!        ONE OR BOTH OF THE ENDPOINTS XLEFT OR XRIGHT, THEN DBDU(K,K)   
!        AND DBDUX(K,K) SHOULD BOTH BE SET TO ZERO WHEN BNDRY IS        
!        CALLED FOR THAT POINT.  WE REFER TO THIS AS A NULL BOUNDARY    
!        CONDITION.  THIS ROUTINE CAN BE STRUCTURED AS FOLLOWS...       
!        THE COMMON BLOCK /ENDPT/ IS NOT A PART OF PDECOL AND           
!        MUST BE SUPPLIED AND DEFINED BY THE USER.                      
!     COMMON /ENDPT/ XLEFT                                              
!     IF( X .NE. XLEFT ) GO TO 10                                       
!        HERE DEFINE AND SET PROPER VALUES FOR DBDU(K,J), DBDUX(K,J),   
!        AND DZDT(K) FOR K,J = 1 TO NPDE FOR THE LEFT BOUNDARY POINT    
!        X = XLEFT.                                                     
!     RETURN                                                            
!  10 CONTINUE                                                          
!        HERE DEFINE AND SET PROPER VALUES FOR DBDU(K,J), DBDUX(K,J),   
!        AND DZDT(K) FOR K,J = 1 TO NPDE FOR THE RIGHT BOUNDARY POINT   
!        X = XRIGHT.                                                    
!     RETURN                                                            
!     END                                                               
!                                                                       
! 3)  SUBROUTINE UINIT( X, U, NPDE )                                    
!     DIMENSION U(NPDE)                                                 
!        THIS ROUTINE IS USED TO PROVIDE THE PDE PACKAGE WITH THE       
!        NEEDED INITIAL CONDITION FUNCTION VALUES.  THE PACKAGE         
!        PROVIDES A VALUE OF THE INPUT VARIABLE X, AND THE USER IS TO   
!        DEFINE THE PROPER INITIAL VALUES (AT T = T0) FOR ALL OF THE    
!        PDE COMPONENTS, I.E.                                           
!           U(K) = DESIRED INITIAL VALUE OF PDE COMPONENT U(K) AT       
!                  X AND T = T0 FOR K = 1 TO NPDE.                      
!        NOTE... THE INCOMING VALUE OF X WILL BE A COLLOCATION POINT    
!        VALUE.  THE INITIAL CONDITIONS AND BOUNDARY CONDITIONS         
!        MUST BE CONSISTENT (SEE ABOVE).                                
!     RETURN                                                            
!     END                                                               
!-----------------------------------------------------------------------
!                                                                       
! OPTIONAL USER SUPPLIED SUBROUTINE                                     
!                                                                       
! IF THE USER DESIRES TO USE THE MF = 11 OR 21 OPTION IN ORDER TO SAVE  
! ABOUT 10-20 PERCENT IN EXECUTION TIME (SEE BELOW), THEN THE USER MUST 
! PROVIDE THE FOLLOWING SUBROUTINE WHICH PROVIDES INFORMATION ABOUT THE 
! DERIVATIVES OF THE FUNCTION F ABOVE. THIS PROVIDES FOR MORE EFFICIENT 
! JACOBIAN MATRIX GENERATION.  ON MOST COMPUTER SYSTEMS, THE USER WILL  
! BE REQUIRED TO SUPPLY THIS SUBROUTINE AS A DUMMY SUBROUTINE IF THE    
! OPTIONS MF = 12 OR 22 ARE USED (SEE BELOW).                           
!                                                                       
! 1)  SUBROUTINE DERIVF( T, X, U, UX, UXX, DFDU, DFDUX, DFDUXX, NPDE )  
!     DIMENSION U(NPDE), UX(NPDE), UXX(NPDE)                            
!     DIMENSION DFDU(NPDE,NPDE), DFDUX(NPDE,NPDE), DFDUXX(NPDE,NPDE)    
!        THE PACKAGE PROVIDES VALUES OF THE INPUT VARIABLES T, X, U, UX,
!        AND UXX, AND THE USER SHOULD CONSTRUCT THIS ROUTINE TO PROVIDE 
!        THE FOLLOWING CORRESPONDING VALUES OF THE OUTPUT ARRAYS        
!        DFDU, DFDUX, AND DFDUXX FOR K,J = 1 TO NPDE...                 
!           DFDU(K,J) = PARTIAL DERIVATIVE OF THE K-TH COMPONENT OF THE 
!                       PDE DEFINING FUNCTION  F  WITH RESPECT TO THE   
!                       VARIABLE U(J).                                  
!          DFDUX(K,J) = PARTIAL DERIVATIVE OF THE K-TH COMPONENT OF THE 
!                       PDE DEFINING FUNCTION  F  WITH RESPECT TO THE   
!                       VARIABLE UX(J).                                 
!         DFDUXX(K,J) = PARTIAL DERIVATIVE OF THE K-TH COMPONENT OF THE 
!                       PDE DEFINING FUNCTION  F  WITH RESPECT TO THE   
!                       VARIABLE UXX(J).                                
!        NOTE... THE INCOMING VALUE OF  X  WILL BE A COLLOCATION POINT  
!        VALUE.                                                         
!     RETURN                                                            
!     END                                                               
!-----------------------------------------------------------------------
!                                                                       
! METHODS USED                                                          
!                                                                       
!   THE PACKAGE PDECOL IS BASED ON THE METHOD OF LINES AND USES A       
!   FINITE ELEMENT COLLOCATION PROCEDURE (WITH PIECEWISE POLYNOMIALS    
!   AS THE TRIAL SPACE) FOR THE DISCRETIZATION OF THE SPATIAL VARIABLE  
!   X.  THE COLLOCATION PROCEDURE REDUCES THE PDE SYSTEM TO A SEMI-     
!   DISCRETE SYSTEM WHICH THEN DEPENDS ONLY ON THE TIME VARIABLE T.     
!   THE TIME INTEGRATION IS THEN ACCOMPLISHED BY USE OF SLIGHTLY        
!   MODIFIED STANDARD TECHNIQUES (SEE REFS. 1,2).                       
!                                                                       
!   PIECEWISE POLYNOMIALS                                               
!                                                                       
!   THE USER IS REQUIRED TO SELECT THE PIECEWISE POLYNOMIAL SPACE       
!   WHICH IS TO BE USED TO COMPUTE HIS APPROXIMATE SOLUTION. FIRST, THE 
!   ORDER, KORD, OF THE POLYNOMIALS TO BE USED MUST BE SPECIFIED        
!   (KORD = POLYNOMIAL DEGREE + 1).  NEXT, THE NUMBER OF PIECES         
!   (INTERVALS), NINT, INTO WHICH THE SPATIAL DOMAIN (XLEFT,XRIGHT) IS  
!   TO BE DIVIDED, IS CHOSEN.  THE NINT + 1 DISTINCT BREAKPOINTS OF     
!   THE DOMAIN MUST BE DEFINED AND SET INTO THE ARRAY XBKPT IN          
!   STRICTLY INCREASING ORDER, I.E.                                     
!   XLEFT=XBKPT(1) .LT. XBKPT(2) .LT. ... .LT. XBKPT(NINT+1)=XRIGHT.    
!   THE APPROXIMATE SOLUTION AT ANY TIME T WILL BE A POLYNOMIAL OF      
!   ORDER KORD OVER EACH SUBINTERVAL (XBKPT(I),XBKPT(I+1)).  THE        
!   NUMBER OF CONTINUITY CONDITIONS, NCC, TO BE IMPOSED ACROSS ALL OF   
!   THE BREAKPOINTS IS THE LAST PIECE OF USER SUPPLIED DATA WHICH IS    
!   REQUIRED TO UNIQUELY DETERMINE THE DESIRED PIECEWISE POLYNOMIAL     
!   SPACE.  FOR EXAMPLE, NCC = 2 WOULD REQUIRE THAT THE APPROXIMATE     
!   SOLUTION (MADE UP OF THE SEPARATE POLYNOMIAL PIECES) AND ITS FIRST  
!   SPATIAL DERIVATIVE BE CONTINUOUS AT THE BREAKPOINTS AND HENCE ON    
!   THE ENTIRE DOMAIN (XLEFT,XRIGHT).  NCC = 3 WOULD REQUIRE THAT THE   
!   APPROXIMATE SOLUTION AND ITS FIRST AND SECOND SPATIAL DERIVATIVES   
!   BE CONTINUOUS AT THE BREAKPOINTS, ETC. THE DIMENSION OF THIS LINEAR 
!   SPACE IS KNOWN AND FINITE AND IS NCPTS = KORD*NINT - NCC*(NINT-1).  
!   THE WELL-KNOWN B-SPLINE BASIS (SEE REF. 3) FOR THIS SPACE IS USED   
!   BY PDECOL AND IT CONSISTS OF NCPTS KNOWN PIECEWISE POLYNOMIAL       
!   FUNCTIONS BF(I,X), FOR I=1 TO NCPTS, WHICH DO NOT DEPEND ON THE     
!   TIME VARIABLE T. WE WISH TO EMPHASIZE THAT THE PIECEWISE POLYNOMIAL 
!   SPACE USED IN PDECOL (WHICH IS SELECTED BY THE USER) WILL DETERMINE 
!   THE MAGNITUDE OF THE SPATIAL DISCRETIZATION ERRORS IN THE COMPUTED  
!   APPROXIMATE SOLUTION.  THE PACKAGE HAS NO CONTROL OVER ERRORS       
!   INTRODUCED BY THE USERS CHOICE OF THIS SPACE.  SEE INPUT PARAMETERS 
!   BELOW.                                                              
!                                                                       
!   COLLOCATION OVER PIECEWISE POLYNOMIALS                              
!                                                                       
!   THE BASIC ASSUMPTION MADE IS THAT THE APPROXIMATE SOLUTION          
!   SATISFIES                                                           
!                       NCPTS                                           
!             U(T,X)  =  SUM  C(I,T) * BF(I,X)                          
!                        I=1                                            
!                                                                       
!   WHERE THE UNKNOWN COEFFICIENTS C DEPEND ONLY ON THE TIME T AND      
!   THE KNOWN BASIS FUNCTIONS DEPEND ONLY ON X (WE HAVE ASSUMED THAT    
!   NPDE = 1 FOR CONVENIENCE).  SO, AT ANY GIVEN TIME T THE APPROX-     
!   IMATE SOLUTION IS A PIECEWISE POLYNOMIAL IN THE USER CHOSEN SPACE.  
!   THE SEMI-DISCRETE EQUATIONS (ACTUALLY ORDINARY DIFFERENTIAL         
!   EQUATIONS) WHICH DETERMINE THE COEFFICIENTS C ARE OBTAINED BY       
!   REQUIRING THAT THE ABOVE APPROXIMATE U(T,X) SATISFY THE PDE AND     
!   BOUNDARY CONDITIONS EXACTLY AT A SET OF NCPTS COLLOCATION POINTS    
!   (SEE COLPNT).  THUS, PDECOL ACTUALLY COMPUTES THE BASIS FUNCTION    
!   COEFFICIENTS RATHER THAN SPECIFIC APPROXIMATE SOLUTION VALUES.      
!                                                                       
!   REFERENCES                                                          
!                                                                       
!   1. MADSEN, N.K. AND R.F. SINCOVEC, PDECOL - COLLOCATION SOFTWARE    
!        FOR PARTIAL DIFFERENTIAL EQUATIONS, ACM-TOMS, VOL.  , NO.  ,   
!                                                                       
!   2. SINCOVEC, R.F. AND N.K. MADSEN, SOFTWARE FOR NONLINEAR PARTIAL   
!        DIFFERENTIAL EQUATIONS, ACM-TOMS, VOL. 1, NO. 3,               
!        SEPTEMBER 1975, PP. 232-260.                                   
!   3. HINDMARSH, A.C., PRELIMINARY DOCUMENTATION OF GEARIB.. SOLUTION  
!        OF IMPLICIT SYSTEMS OF ORDINARY DIFFERENTIAL EQUATIONS WITH    
!        BANDED JACOBIANS, LAWRENCE LIVERMORE LAB, UCID-30130, FEBRUARY 
!        1976.                                                          
!   4. DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.      
!        NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.          
!-----------------------------------------------------------------------
!                                                                       
! USE OF PDECOL                                                         
!                                                                       
! PDECOL IS CALLED ONCE FOR EACH DESIRED OUTPUT VALUE (TOUT) OF THE     
! TIME T, AND IT IN TURN MAKES REPEATED CALLS TO THE CORE INTEGRATOR,   
! STIFIB, WHICH ADVANCES THE TIME BY TAKING SINGLE STEPS UNTIL          
! T .GE. TOUT.  INTERPOLATION TO THE EXACT TIME TOUT IS THEN DONE.      
! SEE TOUT BELOW.                                                       
!                                                                       
!                                                                       
! SUMMARY OF SUGGESTED INPUT VALUES                                     
!                                                                       
!   IT IS OF COURSE IMPOSSIBLE TO SUGGEST INPUT PARAMETER VALUES WHICH  
!   ARE APPROPRIATE FOR ALL PROBLEMS.  THE FOLLOWING SUGGESTIONS ARE TO 
!   BE USED ONLY IF YOU HAVE NO IDEA OF BETTER VALUES FOR YOUR PROBLEM. 
!                                                                       
!   DT    = 1.E-10                                                      
!   XBKPT = CHOOSE  NINT+1  EQUALLY SPACED VALUES SUCH THAT XBKPT(1) =  
!           XLEFT AND  XBKPT(NINT+1) = XRIGHT.                          
!   EPS   = 1.E-4                                                       
!   NINT  = ENOUGH SO THAT ANY FINE STRUCTURE OF THE PROBLEM MAY BE     
!           RESOLVED.                                                   
!   KORD  = 4                                                           
!   NCC   = 2                                                           
!   MF    = 22                                                          
!   INDEX = 1 (ON FIRST CALL ONLY, THEN 0 THEREAFTER).                  
!                                                                       
!                                                                       
! THE INPUT PARAMETERS ARE..                                            
!   T0    =  THE INITIAL VALUE OF T, THE INDEPENDENT VARIABLE           
!              (USED ONLY ON FIRST CALL).                               
!   TOUT  =  THE VALUE OF T AT WHICH OUTPUT IS DESIRED NEXT.  SINCE     
!              THE PACKAGE CHOOSES ITS OWN TIME STEP SIZES, THE         
!              INTEGRATION WILL NORMALLY GO SLIGHTLY BEYOND TOUT        
!              AND THE PACKAGE WILL INTERPOLATE TO T = TOUT.            
!   DT    =  THE INITIAL STEP SIZE IN T, IF INDEX = 1, OR, THE          
!              MAXIMUM STEP SIZE ALLOWED (MUST BE .GT. 0), IF INDEX = 3.
!              USED FOR INPUT ONLY WHEN INDEX = 1 OR 3. SEE BELOW.      
!   XBKPT =  THE ARRAY OF PIECEWISE POLYNOMIAL BREAKPOINTS.             
!              THE NINT+1 VALUES MUST BE STRICTLY INCREASING WITH       
!              XBKPT(1) = XLEFT AND XBKPT(NINT+1) = XRIGHT (USED ONLY   
!              ON FIRST CALL).                                          
!   EPS   =  THE RELATIVE TIME ERROR BOUND  (USED ONLY ON THE           
!              FIRST CALL, UNLESS INDEX = 4).  SINGLE STEP ERROR        
!              ESTIMATES DIVIDED BY CMAX(I) WILL BE KEPT LESS THAN      
!              EPS IN ROOT-MEAN-SQUARE NORM. THE VECTOR CMAX OF WEIGHTS 
!              IS COMPUTED IN PDECOL.  INITIALLY CMAX(I) IS SET TO      
!              ABS(C(I)), WITH A DEFAULT VALUE OF 1 IF ABS(C(I)) .LT. 1.
!              THEREAFTER, CMAX(I) IS THE LARGEST VALUE                 
!              OF ABS(C(I)) SEEN SO FAR, OR THE INITIAL CMAX(I) IF      
!              THAT IS LARGER.  TO ALTER EITHER OF THESE, CHANGE THE    
!              APPROPRIATE STATEMENTS IN THE DO-LOOPS ENDING AT         
!              STATEMENTS 50 AND 130 BELOW.  THE USER SHOULD EXERCISE   
!              SOME DISCRETION IN CHOOSING EPS.  IN GENERAL, THE        
!              OVERALL RUNNING TIME FOR A PROBLEM WILL BE GREATER IF    
!              EPS IS CHOSEN SMALLER. THERE IS USUALLY LITTLE REASON TO 
!              CHOOSE EPS MUCH SMALLER THAN THE ERRORS WHICH ARE BEING  
!              INTRODUCED BY THE USERS CHOICE OF THE POLYNOMIAL SPACE.  
!              SEE RELATED COMMENTS CONCERNING CMAX BELOW STATEMENT 40. 
!   NINT  =  THE NUMBER OF SUBINTERVALS INTO WHICH THE SPATIAL DOMAIN   
!              (XLEFT,XRIGHT) IS TO BE DIVIDED (MUST BE .GE. 1)         
!              (USED ONLY ON FIRST CALL).                               
!   KORD  =  THE ORDER OF THE PIECEWISE POLYNOMIAL SPACE TO BE USED.    
!              ITS VALUE MUST BE GREATER THAN 2 AND LESS THAN 21.  FOR  
!              FIRST ATTEMPTS WE RECOMMEND KORD = 4.  IF THE SOLUTION   
!              IS SMOOTH AND MUCH ACCURACY IS DESIRED, HIGHER VALUES    
!              MAY PROVE TO BE MORE EFFICIENT.  WE HAVE SELDOM USED     
!              VALUES OF KORD IN EXCESS OF 8 OR 9, THOUGH THEY ARE      
!              AVAILABLE FOR USE IN PDECOL (USED ONLY ON FIRST CALL).   
!   NCC   =  THE NUMBER OF CONTINUITY CONDITIONS TO BE IMPOSED ON THE   
!              APPROXIMATE SOLUTION AT THE BREAKPOINTS IN XBKPT.        
!              NCC MUST BE GREATER THAN 1 AND LESS THAN KORD.  WE       
!              RECOMMEND THE USE OF NCC = 2 (WITH NOGAUS = 0, SEE       
!              BELOW), SINCE THEORY PREDICTS THAT DRAMATICALLY MORE     
!              ACCURATE RESULTS CAN OFTEN BE OBTAINED USING THIS CHOICE 
!              (USED ONLY ON FIRST CALL).                               
!   NPDE  =  THE NUMBER OF PARTIAL DIFFERENTIAL EQUATIONS IN THE SYSTEM 
!              TO BE SOLVED (USED ONLY ON FIRST CALL).                  
!   MF    =  THE METHOD FLAG  (USED ONLY ON FIRST CALL, UNLESS          
!              INDEX = 4).  ALLOWED VALUES ARE 11, 12, 21, 22.          
!              FOR FIRST ATTEMPTS WE RECOMMEND THE USE OF MF = 22.      
!              MF HAS TWO DECIMAL DIGITS, METH AND MITER                
!              (MF = 10*METH + MITER).                                  
!              METH IS THE BASIC METHOD INDICATOR..                     
!                METH = 1  MEANS THE ADAMS METHODS (GENERALIZATIONS OF  
!                          CRANK-NICOLSON).                             
!                METH = 2  MEANS THE BACKWARD DIFFERENTIATION           
!                          FORMULAS (BDF), OR STIFF METHODS OF GEAR.    
!              MITER IS THE ITERATION METHOD INDICATOR                  
!              AND DETERMINES HOW THE JACOBIAN MATRIX IS                
!              TO BE COMPUTED..                                         
!                MITER = 1 MEANS CHORD METHOD WITH ANALYTIC JACOBIAN.   
!                          FOR THIS USER SUPPLIES SUBROUTINE DERIVF.    
!                          SEE DESCRIPTION ABOVE.                       
!                MITER = 2 MEANS CHORD METHOD WITH JACOBIAN CALCULATED  
!                          INTERNALLY BY FINITE DIFFERENCES.  SEE       
!                          SUBROUTINES PSETIB AND DIFFF.                
!   INDEX =  INTEGER USED ON INPUT TO INDICATE TYPE OF CALL,            
!              WITH THE FOLLOWING VALUES AND MEANINGS..                 
!                 1    THIS IS THE FIRST CALL FOR THIS PROBLEM.         
!                 0    THIS IS NOT THE FIRST CALL FOR THIS PROBLEM,     
!                      AND INTEGRATION IS TO CONTINUE.                  
!                 2    SAME AS 0 EXCEPT THAT TOUT IS TO BE HIT          
!                      EXACTLY (NO INTERPOLATION IS DONE).  SEE NOTE    
!                      BELOW.  ASSUMES TOUT .GE. THE CURRENT T.         
!                      IF TOUT IS .LT. THE CURRENT TIME, THEN TOUT IS   
!                      RESET TO THE CURRENT TIME AND CONTROL IS         
!                      RETURNED TO THE USER.  A CALL TO VALUES WILL     
!                      PRODUCE ANSWERS FOR THE NEW VALUE OF TOUT.       
!                 3    SAME AS 0 EXCEPT CONTROL RETURNS TO CALLING      
!                      PROGRAM AFTER ONE STEP.  TOUT IS IGNORED AND     
!                      DT MUST BE SET .GT. 0 TO A MAXIMUM ALLOWED       
!                      DT VALUE. SEE ABOVE.                             
!                 4    THIS IS NOT THE FIRST CALL FOR THE PROBLEM,      
!                      AND THE USER HAS RESET EPS AND/OR MF.            
!              SINCE THE NORMAL OUTPUT VALUE OF INDEX IS 0,             
!              IT NEED NOT BE RESET FOR NORMAL CONTINUATION.            
!                                                                       
! NOTE.. THE PACKAGE MUST HAVE TAKEN AT LEAST ONE SUCCESSFUL TIME       
! STEP BEFORE A CALL WITH INDEX = 2 OR 4 IS ALLOWED.                    
! AFTER THE INITIAL CALL, IF A NORMAL RETURN OCCURRED AND A NORMAL      
! CONTINUATION IS DESIRED, SIMPLY RESET TOUT AND CALL AGAIN.            
! ALL OTHER PARAMETERS WILL BE READY FOR THE NEXT CALL.                 
! A CHANGE OF PARAMETERS WITH INDEX = 4 CAN BE MADE AFTER               
! EITHER A SUCCESSFUL OR AN UNSUCCESSFUL RETURN PROVIDED AT LEAST       
! ONE SUCCESSFUL TIME STEP HAS BEEN MADE.                               
!                                                                       
!   WORK  =  FLOATING POINT WORKING ARRAY FOR PDECOL.  WE RECOMMEND     
!              THAT IT BE INITIALIZED TO ZERO BEFORE THE FIRST CALL     
!              TO PDECOL.  ITS TOTAL LENGTH MUST BE AT LEAST            
!                                                                       
!              KORD + 4*NPDE + 9*NPDE**2 + NCPTS*(3*KORD + 2) +         
!              NPDE*NCPTS*(3*ML + MAXDER + 7)                           
!                                                                       
!              WHERE ML AND MAXDER ARE DEFINED BELOW (SEE STORAGE       
!              ALLOCATION).                                             
!   IWORK =  INTEGER WORKING ARRAY FOR PDECOL.  THE FIRST TWO           
!              LOCATIONS MUST BE DEFINED AS FOLLOWS...                  
!              IWORK(1) = LENGTH OF USERS ARRAY WORK                    
!              IWORK(2) = LENGTH OF USERS ARRAY IWORK                   
!              THE TOTAL LENGTH OF IWORK MUST BE AT LEAST               
!              NCPTS*(NPDE + 1).                                        
! OUTPUT                                                                
!                                                                       
! THE SOLUTION VALUES ARE NOT RETURNED DIRECTLY TO THE USER BY PDECOL.  
! THE METHODS USED IN PDECOL COMPUTE BASIS FUNCTION COEFFICIENTS, SO    
! THE USER (AFTER A RETURN FROM PDECOL) MUST CALL THE PACKAGE ROUTINE   
! VALUES TO OBTAIN HIS APPROXIMATE SOLUTION VALUES AT ANY DESIRED SPACE 
! POINTS X AT THE TIME T = TOUT.  SEE THE COMMENTS IN SUBROUTINE VALUES 
! FOR DETAILS ON HOW TO PROPERLY MAKE THE CALL.                         
!                                                                       
! THE COMMON BLOCK /GEAR0/ CAN BE ACCESSED EXTERNALLY BY THE USER       
! IF DESIRED.  IT CONTAINS THE STEP SIZE LAST USED (SUCCESSFULLY),      
! THE ORDER LAST USED (SUCCESSFULLY), THE NUMBER OF STEPS TAKEN         
! SO FAR, THE NUMBER OF RESIDUAL EVALUATIONS (RES CALLS) SO FAR,        
! AND THE NUMBER OF MATRIX EVALUATIONS (PSETIB CALLS) SO FAR.           
! DIFFUN CALLS ARE COUNTED IN WITH RESIDUAL EVALUATIONS.                
!                                                                       
! THE OUTPUT PARAMETERS ARE..                                           
!   DT    =  THE STEP SIZE USED LAST, WHETHER SUCCESSFULLY OR NOT.      
!   TOUT  =  THE OUTPUT VALUE OF T.  IF INTEGRATION WAS SUCCESSFUL,     
!              AND THE INPUT VALUE OF INDEX WAS NOT 3, TOUT IS          
!              UNCHANGED FROM ITS INPUT VALUE.  OTHERWISE, TOUT         
!              IS THE CURRENT VALUE OF T TO WHICH THE INTEGRATION       
!              HAS BEEN COMPLETED.                                      
!   INDEX =  INTEGER USED ON OUTPUT TO INDICATE RESULTS,                
!              WITH THE FOLLOWING VALUES AND MEANINGS..                 
!         0    INTEGRATION WAS COMPLETED TO TOUT OR BEYOND.             
!        -1    THE INTEGRATION WAS HALTED AFTER FAILING TO PASS THE     
!              ERROR TEST EVEN AFTER REDUCING DT BY A FACTOR OF         
!              1.E10 FROM ITS INITIAL VALUE.                            
!        -2    AFTER SOME INITIAL SUCCESS, THE INTEGRATION WAS          
!              HALTED EITHER BY REPEATED ERROR TEST FAILURES OR BY      
!              A TEST ON EPS.  TOO MUCH ACCURACY HAS BEEN REQUESTED.    
!        -3    THE INTEGRATION WAS HALTED AFTER FAILING TO ACHIEVE      
!              CORRECTOR CONVERGENCE EVEN AFTER REDUCING DT BY A        
!              FACTOR OF 1.E10 FROM ITS INITIAL VALUE.                  
!        -4    SINGULAR MATRIX ENCOUNTERED.  PROBABLY DUE TO STORAGE    
!              OVERWRITES.                                              
!        -5    INDEX WAS 4 ON INPUT, BUT THE DESIRED CHANGES OF         
!              PARAMETERS WERE NOT IMPLEMENTED BECAUSE TOUT             
!              WAS NOT BEYOND T.  INTERPOLATION TO T = TOUT WAS         
!              PERFORMED AS ON A NORMAL RETURN.  TO TRY AGAIN,          
!              SIMPLY CALL AGAIN WITH INDEX = 4 AND A NEW TOUT.         
!        -6    ILLEGAL INDEX VALUE.                                     
!        -7    ILLEGAL EPS VALUE.                                       
!        -8    AN ATTEMPT TO INTEGRATE IN THE WRONG DIRECTION.  THE     
!              SIGN OF DT IS WRONG RELATIVE TO T0 AND TOUT.             
!        -9    DT .EQ. 0.0.                                             
!       -10    ILLEGAL NINT VALUE.                                      
!       -11    ILLEGAL KORD VALUE.                                      
!       -12    ILLEGAL NCC VALUE.                                       
!       -13    ILLEGAL NPDE VALUE.                                      
!       -14    ILLEGAL MF VALUE.                                        
!       -15    ILLEGAL BREAKPOINTS - NOT STRICTLY INCREASING.           
!       -16    INSUFFICIENT STORAGE FOR WORK OR IWORK.                  
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
! SUMMARY OF ALL PACKAGE ROUTINES                                       
!                                                                       
! PDECOL - STORAGE ALLOCATION, ERROR CHECKING, INITIALIZATION, REPEATED 
!          CALLS TO STIFIB TO ADVANCE THE TIME.                         
!                                                                       
! INTERP - INTERPOLATES COMPUTED BASIS FUNCTION COEFFICIENTS TO THE     
!          DESIRED OUTPUT TIMES, TOUT, FOR USE BY VALUES.               
!                                                                       
! INITAL - INITIALIZATION, GENERATION AND STORAGE OF PIECEWISE POLY-    
!          NOMIAL SPACE BASIS FUNCTION VALUES AND DERIVATIVES, DET-     
!          ERMINES THE BASIS FUNCTION COEFFICINTS OF THE PIECEWISE      
!          POLYNOMIALS WHICH INTERPOLATE THE USERS INITIAL CONDITIONS.  
!                                                                       
! COLPNT - GENERATION OF REQUIRED COLLOCATION POINTS.                   
!                                                                       
! BSPLVD - B-SPLINE PACKAGE ROUTINES WHICH ALLOW FOR EVALUATION OF      
! BSPLVN   ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.             
! INTERV                                                                
!                                                                       
! VALUES - GENERATION AT ANY POINT(S) OF VALUES OF THE COMPUTED         
!          APPROXIMATE SOLUTION AND ITS DERIVATIVES WHICH ARE           
!          PIECEWISE POLYNOMIALS.  THE SUBROUTINE IS CALLED ONLY BY     
!          THE USER.                                                    
!                                                                       
! STIFIB - CORE INTEGRATOR, TAKES SINGLE TIME STEPS TO ADVANCE THE      
!          TIME.  ASSEMBLES AND SOLVES THE PROPER NONLINEAR EQUATIONS   
!          WHICH ARE RELATED TO USE OF ADAMS OR GEAR TYPE INTEGRATION   
!          FORMULAS.  CHOOSES PROPER STEP SIZE AND INTEGRATION FORMULA  
!          ORDER TO MAINTAIN A DESIRED ACCURACY.  DESIGNED FOR ODE      
!          PROBLEMS OF THE FORM A * (DY/DT) = G(T,Y).                   
!                                                                       
! CSTCOL  - GENERATES INTEGRATION FORMULA AND ERROR CONTROL COEFFICIENTS
!                                                                       
! RES    - COMPUTES RESIDUAL VECTORS USED IN SOLVING THE NONLINEAR      
!          EQUATIONS BY A MODIFIED NEWTON METHOD.                       
!                                                                       
! DIFFUN - COMPUTES A**-1 * G(T,Y) WHERE A AND G ARE AS ABOVE (STIFIB). 
!                                                                       
! ADDA   - ADDS THE A MATRIX TO A GIVEN MATRIX IN BAND FORM.            
!                                                                       
! EVAL   - EVALUATES THE COMPUTED PIECEWISE POLYNOMIAL SOLUTION AND     
!          DERIVATIVES AT COLLOCATION POINTS.                           
!                                                                       
! GFUN   - EVALUATES THE FUNCTION G(T,Y) BY CALLING EVAL AND THE USER   
!          SUBROUTINES F AND BNDRY.                                     
!                                                                       
! PSETIB - GENERATES PROPER JACOBIAN MATRICES REQUIRED BY THE MODIFIED  
!          NEWTON METHOD.                                               
!                                                                       
! DIFFF  - PERFORMS SAME ROLE AS THE USER ROUTINE DERIVF. COMPUTES      
!          DERIVATIVE APPROXIMATIONS BY USE OF FINITE DIFFERENCES.      
!                                                                       
! DECB   - PERFORM AN LU DECOMPOSTION AND FORWARD AND BACKWARD          
! SOLB     SUBSTITUTION FOR SOLVING BANDED SYSTEMS OF LINEAR EQUATIONS. 
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
! STORAGE ALLOCATION                                                    
!                                                                       
! SINCE PDECOL IS A DYNAMICALLY DIMENSIONED PROGRAM, MOST OF ITS        
! WORKING STORAGE IS PROVIDED BY THE USER IN THE ARRAYS WORK AND IWORK. 
! THE FOLLOWING GIVES A LIST OF THE ARRAYS WHICH MAKE UP THE CONTENTS   
! WORK AND IWORK, THEIR LENGTHS, AND THEIR USES.  WHEN MORE THAN ONE    
! NAME IS GIVEN, IT INDICATES THAT DIFFERENT NAMES ARE USED FOR THE     
! SAME ARRAY IN DIFFERENT PARTS OF THE PROGRAM.  THE DIFFERENT NAMES    
! OCCUR BECAUSE PDECOL IS AN AMALGAMATION OF SEVERAL OTHER CODES        
! WRITTEN BY DIFFERENT PEOPLE AND WE HAVE TRIED TO LEAVE THE SEPARATE   
! PARTS AS UNCHANGED FROM THEIR ORIGINAL VERSIONS AS POSSIBLE.          
!                                                                       
!                                                                       
!   NAMES       LENGTH        USE                                       
! ---------     ------------  -------------------------------------     
!                                                                       
! BC            4*NPDE**2     BOUNDARY CONDITION INFORMATION.           
! WORK                                                                  
!                                                                       
! A             3*KORD*NCPTS  BASIS FUNCTION VALUES AT COLLOCATION POINT
! WORK(IW1)                                                             
!                                                                       
! XT            NCPTS + KORD  BREAKPOINT SEQUENCE FOR GENERATION OF BASI
! WORK(IW2)                   FUNCTION VALUES.                          
!                                                                       
! XC            NCPTS         COLLOCATION POINTS.                       
! WORK(IW3)                                                             
!                                                                       
! CMAX          NPDE*NCPTS    VALUES USED IN ESTIMATING TIME            
! YMAX                        INTEGRATION ERRORS.                       
! WORK(IW4)                                                             
!                                                                       
! ERROR         NPDE*NCPTS    TIME INTEGRATION ERRORS.                  
! WORK(IW5)                                                             
!                                                                       
! SAVE1         NPDE*NCPTS    WORKING STORAGE FOR THE TIME INTEGRATION  
! WORK(IW6)                   METHOD.                                   
!                                                                       
! SAVE2         NPDE*NCPTS    WORKING STORAGE FOR THE TIME INTEGRATION  
! WORK(IW7)                   METHOD.                                   
!                                                                       
! SAVE3         NPDE*NCPTS    WORKING STORAGE FOR THE TIME INTEGRATION  
! WORK(IW8)                   METHOD.                                   
!                                                                       
! UVAL          3*NPDE        WORKING STORAGE FOR VALUES OF U, UX, AND  
! WORK(IW9)                   USS AT ONE POINT.                         
!                                                                       
! C             NPDE*NCPTS*   CURRENT BASIS FUNCTION COEFFICIENT VALUES 
! Y               (MAXDER+1)  AND THEIR SCALED TIME DERIVATIVES.        
! WORK(IW10)                                                            
!                                                                       
! DFDU          NPDE**2       WORKING STORAGE USED TO COMPUTE THE       
! WORK(IW11)                  JACOBIAN MATRIX.                          
!                                                                       
! DFDUX         NPDE**2       WORKING STORAGE USED TO COMPUTE THE       
! WORK(IW12)                  JACOBIAN MATRIX.                          
!                                                                       
! DFDUXX        NPDE**2       WORKING STORAGE USED TO COMPUTE THE       
! WORK(IW13)                  JACOBIAN MATRIX.                          
!                                                                       
! DBDU          NPDE**2       BOUNDARY CONDITION INFORMATION.           
! WORK(IW14)                                                            
!                                                                       
! DBDUX         NPDE**2       BOUNDARY CONDITION INFORMATION.           
! WORK(IW15)                                                            
!                                                                       
! DZDT          NPDE          BOUNDARY CONDITION INFORMATION.           
! WORK(IW16)                                                            
!                                                                       
! PW            NPDE*NCPTS*   STORAGE AND PROCESSING OF THE JACOBIAN    
! WORK(IW17)      (3*ML+1)    MATRIX.                                   
!                                                                       
! ILEFT         NCPTS         POINTERS TO BREAKPOINT SEQUENCE FOR       
! IWORK                       GENERATION OF BASIS FUNCTION VALUES.      
!                                                                       
! IPIV          NPDE*NCPTS    PIVOT INFORMATION FOR THE LU DECOMPOSED   
! IWORK(IW18)                 JACOBIAN MATRIX PW.                       
!                                                                       
! WHERE...                                                              
!                                                                       
!      NCPTS = KORD*NINT - NCC*(NINT-1)                                 
!         ML = NPDE*(KORD+IQUAD-1) - 1                                  
!      IQUAD = 1 IF KORD = 3 AND A NULL BOUNDARY CONDITION EXISTS       
!      IQUAD = 0 OTHERWISE                                              
!     MAXDER = 5 UNLESS OTHERWISE SET BY THE USER INTO /OPTION/.        
!                                                                       
! THE COMMON BLOCK /OPTION/ CONTAINS THE VARIABLES NOGAUS AND MAXDER.   
! NOGAUS IS SET .EQ. 0 IN THE BLOCK DATA.  IT CAN BE CHANGED TO BE      
! SET .EQ. 1 IF THE GAUSS-LEGENDRE COLLOCATION POINTS ARE NOT           
! DESIRED WHEN NCC = 2 (SEE ABOVE AND COLPNT).  MAXDER IS SET           
! .EQ. 5 IN THE BLOCK DATA AND ITS VALUE REPRESENTS THE                 
! MAXIMUM ORDER OF TIME INTEGRATION FORMULA ALLOWED.  ITS VALUE         
! AFFECTS THE STORAGE REQUIRED IN WORK AND MAY BE CHANGED IF            
! DESIRED.  SEE CSTCOL FOR RESTRICTIONS.  THESE CHANGES MAY BE MADE BY  
! THE USER BY ACCESSING /OPTION/ IN HIS CALLING PROGRAM (BEFORE THE     
! FIRST CALL TO PDECOL) OR BY CHANGING THE DATA STATEMENT IN            
! THE BLOCK DATA.                                                       
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
! COMMUNICATION                                                         
!                                                                       
! EACH SUBROUTINE IN THE PACKAGE CONTAINS A COMMUNICATION SUMMARY       
! AS INDICATED BELOW.                                                   
!                                                                       
! PACKAGE ROUTINES CALLED..   EVAL,INITAL,INTERP,STIFIB                 
! USER ROUTINES CALLED..      BNDRY                                     
! CALLED BY..                 USERS MAIN PROGRAM                        
! FORTRAN FUNCTIONS USED..    ABS,AMAX1,FLOAT,SQRT                      
!-----------------------------------------------------------------------
      DIMENSION WORK ( * ), IWORK ( * ), XBKPT ( * ) 
      COMMON / GEAR0 / DTUSED, NQUSED, NSTEP, NFE, NJE 
      COMMON / GEAR1 / T, DTC, DTMN, DTMX, EPSC, UROUND, N, MFC, KFLAG, &
      JSTART                                                            
      COMMON / GEAR9 / EPSJ, R0, ML, MU, MW, NM1, N0ML, N0W 
!old      COMMON / OPTION / NOGAUS, MAXDER 
      COMMON / SIZES / NIN, KOR, NC, NPD, NCPTS, NEQN, IQUAD 
!old      COMMON / ISTART / IW1, IW2, IW3, IW4, IW5, IW6, IW7, IW8, IW9,    &
!old   IW10, IW11, IW12, IW13, IW14, IW15, IW16, IW17, IW18              
      COMMON / IOUNIT / LOUT 
!new(                                                                   
      INTEGER nDT 
!new)                                                                   
!dbg(
!      real UVAL(NPDE,3)
!dbg)
      IF (INDEX.NE.1) GOTO 5 
!      SET VARIABLES HERE INSTEAD OF BLOCK DATA TO AVOID LOADER PROBLEMS
!    LOUT   = THE LOGICAL UNIT NUMBER FOR THE OUTPUT OF MESSAGES DURING 
!             THE INTEGRATION.                                          
!    NOGAUS = SET .EQ. 1 IF THE GAUSS-LEGENDRE COLLOCATION POINTS ARE   
!             NOT DESIRED WHEN NCC = 2 (SEE PDECOL AND COLPNT).         
!    MAXDER = SET .EQ. 5.  ITS VALUE REPRESENTS THE MAXIMUM ORDER OF    
!             THE TIME INTEGRATION ALLOWED.  ITS VALUE AFFECTS THE STOR-
!             AGE REQUIRED IN /HISTRY/ AND MAY BE CHANGED IF DESIRED    
!             (SEE CSTCOL FOR RESTRICTIONS).                            
      LOUT = 6 
!dbg(                                                                   
!                                                                       
!     The Gauss-Legendre collocation method gives a corrupt solution!!! 
!                                                                       
!???        NOGAUS=0                                                    
!old      NOGAUS = 1 
!dbg)                                                                   
!old      MAXDER = 5 
      UROUND = R1MACH (4) 
    5 CONTINUE 
      IF (INDEX.EQ.0) GOTO 60 
      IF (INDEX.EQ.2) GOTO 70 
      IF (INDEX.EQ.4) GOTO 80 
      IF (INDEX.EQ.3) GOTO 90 
!-----------------------------------------------------------------------
! SEVERAL CHECKS ARE MADE HERE TO DETERMINE IF THE INPUT PARAMETERS     
! HAVE LEGAL VALUES.  ERROR CHECKS ARE MADE ON INDEX, EPS, (T0-TOUT)*DT,
! DT, NINT, KORD, NCC, NPDE, MF, WHETHER THE BREAKPOINT SEQUENCE IS     
! STRICTLY INCREASING, AND WHETHER THERE IS SUFFICIENT STORAGE          
! PROVIDED FOR WORK AND IWORK.  PROBLEM DEPENDENT PARAMETERS ARE        
! CALCULATED AND PLACED IN COMMON.                                      
!-----------------------------------------------------------------------
      IERID = - 6 
      IF (INDEX.NE.1) GOTO 320 
      IERID = IERID-1 
      IF (EPS.LE.0.) GOTO 320 
      IERID = IERID-1 
      IF ( (T0 - TOUT) * DT.GT.0.) GOTO 320 
      IERID = IERID-1 
      IF (DT.EQ.0.0) GOTO 320 
      IERID = IERID-1 
      NIN = NINT 
      IF (NIN.LT.1) GOTO 320 
      IERID = IERID-1 
      KOR = KORD 
      IF (KOR.LT.3.OR.KOR.GT.20) GOTO 320 
      IERID = IERID-1 
      NC = NCC 
      IF (NCC.LT.2.OR.NCC.GE.KOR) GOTO 320 
      IERID = IERID-1 
      NPD = NPDE 
      NPDE2 = NPD * NPD 
      IF (NPDE.LT.1) GOTO 320 
      IERID = IERID-1 
      IF (MF.NE.22.AND.MF.NE.21.AND.MF.NE.12.AND.MF.NE.11) GOTO 320 
      IERID = IERID-1 
      DO 10 K = 1, NIN 
        IF (XBKPT (K) .GE.XBKPT (K + 1) ) GOTO 320 
   10 END DO 
      NCPTS = KOR + (NIN - 1) * (KOR - NCC) 
      NEQN = NPDE * NCPTS 
      ML = (KOR - 1) * NPDE-1 
      MU = ML 
      MW = ML + ML + 1 
      N0W = NEQN * MW 
      IWSAVE = IWORK (1) 
      IISAVE = IWORK (2) 
      IW1 = 4 * NPDE2 + 1 
      IW2 = IW1 + 3 * KORD * NCPTS 
      IW3 = IW2 + NCPTS + KORD 
      IW4 = IW3 + NCPTS 
      IW5 = IW4 + NEQN 
      IW6 = IW5 + NEQN 
      IW7 = IW6 + NEQN 
      IW8 = IW7 + NEQN 
      IW9 = IW8 + NEQN 
      IW10 = IW9 + 3 * NPDE 
      IW11 = IW10 + NEQN * (MAXDER + 1) 
      IW12 = IW11 + NPDE2 
      IW13 = IW12 + NPDE2 
      IW14 = IW13 + NPDE2 
      IW15 = IW14 + NPDE2 
      IW16 = IW15 + NPDE2 
      IW17 = IW16 + NPDE 
      IW18 = NCPTS + 1 
      IERID = IERID-1 
      IWSTOR = IW17 + NEQN * (3 * ML + 1) - 1 
      IISTOR = IW18 + NEQN - 1 
      IF (IWSAVE.LT.IWSTOR.OR.IISAVE.LT.IISTOR) GOTO 335 
!-----------------------------------------------------------------------
! PERFORM INITIALIZATION TASKS.  IF KORD .EQ. 3 THEN CALCULATE THE BAND-
! WIDTH OF THE ASSOCIATED MATRIX PROBLEM BY DETERMINING THE TYPE OF     
! BOUNDARY CONDITIONS, THEN CHECK FOR SUFFICIENT STORAGE AGAIN.         
!-----------------------------------------------------------------------
      CALL INITAL (KOR, WORK (IW1), WORK (IW6), XBKPT, WORK (IW2),      &
      WORK (IW3), WORK (IW17), IWORK (IW18), IWORK)                     
!dbg(
!      do i = 1,NCPTS
!         call EVAL(i,NPDE,WORK(IW6),UVAL,WORK(IW1),IWORK)
!         write (*,'(i3,1p,3(1x,e9.3))'),i,WORK(IW3+i-1),UVAL(1,1),UVAL(2,1)
!      end do
!      write (*,*) "USING INITIAL SOLUTION FROM FILE"
!dbg)
      IF (IQUAD.NE.0) GOTO 280 
      IF (KOR.NE.3) GOTO 40 
      CALL EVAL (1, NPDE, WORK (IW6), WORK (IW9), WORK (IW1), IWORK) 
!old      CALL BNDRY(T0,WORK(IW3),WORK(IW9),WORK(IW9+NPDE),WORK(IW14),  
!old     *           WORK(IW15),WORK(IW16),NPDE)                        
      CALL BNDRY (T0, WORK (IW3), WORK (IW9), WORK (IW9 + NPDE),        &
      WORK (IW14), WORK (IW15), WORK (IW16) )                           
      DO 20 K = 1, NPDE 
        I = K + NPDE * (K - 1) - 1 
        IF (WORK (IW14 + I) .EQ.0.0.AND.WORK (IW15 + I) .EQ.0.0) IQUAD =&
        1                                                               
   20 END DO 
      CALL EVAL (NCPTS, NPDE, WORK (IW6), WORK (IW9), WORK (IW1),       &
      IWORK)                                                            
!old      CALL BNDRY(T0,WORK(IW3+NCPTS-1),WORK(IW9),WORK(IW9+NPDE),     
!old     *           WORK(IW14),WORK(IW15),WORK(IW16),NPDE)             
      CALL BNDRY (T0, WORK (IW3 + NCPTS - 1), WORK (IW9), WORK (IW9 +   &
      NPDE), WORK (IW14), WORK (IW15), WORK (IW16) )                    
      DO 30 K = 1, NPDE 
        I = K + NPDE * (K - 1) - 1 
        IF (WORK (IW14 + I) .EQ.0.0.AND.WORK (IW15 + I) .EQ.0.0) IQUAD =&
        1                                                               
   30 END DO 
      ML = ML + IQUAD * NPDE 
      MU = ML 
      MW = ML + ML + 1 
      N0W = NEQN * MW 
   40 CONTINUE 
      IWSTOR = IW17 + NEQN * (3 * ML + 1) - 1 
      IF (IWSAVE.LT.IWSTOR) GOTO 335 
!-----------------------------------------------------------------------
! IF INITIAL VALUES OF CMAX OTHER THAN THOSE SET BELOW ARE DESIRED,     
! THEY SHOULD BE SET HERE.  ALL CMAX(I) MUST BE POSITIVE.               
! HAVING PROPER VALUES OF CMAX FOR THE PROBLEM BEING SOLVED IS AS       
! IMPORTANT AS CHOOSING EPS (SEE ABOVE), SINCE ERRORS ARE               
! MEASURED RELATIVE TO CMAX.  IF VALUES FOR DTMN OR DTMX, THE           
! BOUNDS ON ABS(DT), OTHER THAN THOSE BELOW ARE DESIRED, THEY           
! SHOULD BE SET BELOW.                                                  
!-----------------------------------------------------------------------
      DO 50 I = 1, NEQN 
        I1 = I - 1 
        WORK (IW4 + I1) = ABS (WORK (IW6 + I1) ) 
        IF (WORK (IW4 + I1) .LT.1.) WORK (IW4 + I1) = 1. 
   50 WORK (IW10 + I1) = WORK (IW6 + I1) 
      N = NEQN 
      T = T0 
      DTC = DT 
      DTMN = ABS (DT) 
      DTUSED = 0. 
      EPSC = EPS 
      MFC = MF 
      JSTART = 0 
      EPSJ = SQRT (UROUND) 
      NM1 = NEQN - 1 
      N0ML = NEQN * ML 
      NHCUT = 0 
      KFLAG = 0 
      TOUTP = T0 
      IF (T0.EQ.TOUT) GOTO 360 
   60 DTMX = ABS (TOUT - TOUTP) * 10. 
      GOTO 140 
!                                                                       
   70 DTMX = ABS (TOUT - TOUTP) * 10. 
      IF ( (T - TOUT) * DTC.GE.0.) GOTO 340 
      GOTO 150 
!                                                                       
   80 IF ( (T - TOUT) * DTC.GE.0.) GOTO 300 
      JSTART = - 1 
      EPSC = EPS 
      MFC = MF 
      GOTO 100 
!                                                                       
   90 DTMX = DT 
!                                                                       
!new(                                                                   
      print *,"KUKU"
      nDT = 0 
  100 CONTINUE 
!new(                                                                   
      IF (mod (nDT, 1000) .eq.0) print '("Time loop T1 = ",1p,e9.3,&
     &" n = ",i7," T = ",e9.3,                " DT = ",e9.3)', TOUT, nDT, T, D&
     &TC                                                                
      nDT = nDT + 1 
!new)                                                                   
      IF ( (T + DTC) .EQ.T) WRITE (LOUT, '(" WARNING..  T + DT = T ON NE&
     &XT STEP.")')                                                      
!old  110 FORMAT(36H WARNING..  T + DT = T ON NEXT STEP.)               
!-----------------------------------------------------------------------
! TAKE A TIME STEP BY CALLING THE INTEGRATOR.                           
!-----------------------------------------------------------------------
      CALL STIFIB (NEQN, WORK (IW10), WORK (IW4), WORK (IW5), WORK (IW6)&
      , WORK (IW7), WORK (IW8), WORK (IW17), IWORK (IW18), WORK, IWORK) 
!                                                                       
      KGO = 1 - KFLAG 
      GOTO (120, 160, 220, 260, 280), KGO 
! KFLAG  =   0,  -1,  -2,  -3   -4                                      
!                                                                       
  120 CONTINUE 
!-----------------------------------------------------------------------
! NORMAL RETURN FROM INTEGRATOR.                                        
!                                                                       
! THE WEIGHTS CMAX(I) ARE UPDATED.  IF DIFFERENT VALUES ARE DESIRED,    
! THEY SHOULD BE SET HERE.  A TEST IS MADE FOR EPS BEING TOO SMALL      
! FOR THE MACHINE PRECISION.                                            
!                                                                       
! ANY OTHER TESTS OR CALCULATIONS THAT ARE REQUIRED AFTER EVERY         
! STEP SHOULD BE INSERTED HERE.                                         
!                                                                       
! IF INDEX = 3, SAVE1 IS SET TO THE CURRENT C VALUES ON RETURN.         
! IF INDEX = 2, DT IS CONTROLLED TO HIT TOUT (WITHIN ROUNDOFF           
! ERROR), AND THEN THE CURRENT C VALUES ARE PUT IN SAVE1 ON RETURN.     
! FOR ANY OTHER VALUE OF INDEX, CONTROL RETURNS TO THE INTEGRATOR       
! UNLESS TOUT HAS BEEN REACHED.  THEN INTERPOLATED VALUES OF C ARE      
! COMPUTED AND STORED IN SAVE1 ON RETURN.                               
! IF INTERPOLATION IS NOT DESIRED, THE CALL TO INTERP SHOULD BE         
! REMOVED AND CONTROL TRANSFERRED TO STATEMENT 340 INSTEAD OF 360.      
!-----------------------------------------------------------------------
      D = 0. 
      DO 130 I = 1, NEQN 
        I1 = I - 1 
        AYI = ABS (WORK (IW10 + I1) ) 
        WORK (IW4 + I1) = AMAX1 (WORK (IW4 + I1), AYI) 
  130 D = D+ (AYI / WORK (IW4 + I1) ) **2 
      D = D * (UROUND / EPS) **2 
      IF (D.GT.FLOAT (NEQN) ) GOTO 240 
      IF (INDEX.EQ.3) GOTO 340 
      IF (INDEX.EQ.2) GOTO 150 
  140 IF ( (T - TOUT) * DTC.LT.0.) GOTO 100 
      CALL INTERP (TOUT, WORK (IW10), NEQN, WORK (IW6) ) 
      GOTO 360 
!                                                                       
  150 IF ( ( (T + DTC) - TOUT) * DTC.LE.0.) GOTO 100 
      IF (ABS (T - TOUT) .LE.100. * UROUND * DTMX) GOTO 340 
      IF ( (T - TOUT) * DTC.GE.0.) GOTO 340 
      DTC = (TOUT - T) * (1. - 4. * UROUND) 
      JSTART = - 1 
      GOTO 100 
!-----------------------------------------------------------------------
! ON AN ERROR RETURN FROM INTEGRATOR, AN IMMEDIATE RETURN OCCURS IF     
! KFLAG = -2 OR -4, AND RECOVERY ATTEMPTS ARE MADE OTHERWISE.           
! TO RECOVER, DT AND DTMN ARE REDUCED BY A FACTOR OF .1 UP TO 10        
! TIMES BEFORE GIVING UP.                                               
!-----------------------------------------------------------------------
  160 WRITE (LOUT, 170) T 
  170 FORMAT(//35H KFLAG = -1 FROM INTEGRATOR AT T = ,E16.8/            &
     &   40H  ERROR TEST FAILED WITH ABS(DT) = DTMIN/)                  
  180 IF (NHCUT.EQ.10) GOTO 200 
      NHCUT = NHCUT + 1 
      DTMN = .1 * DTMN 
      DTC = .1 * DTC 
      WRITE (LOUT, 190) DTC 
  190 FORMAT(25H  DT HAS BEEN REDUCED TO ,E16.8,                        &
     &   26H  AND STEP WILL BE RETRIED//)                               
      JSTART = - 1 
      GOTO 100 
!                                                                       
  200 WRITE (LOUT, 210) 
  210 FORMAT(//44H PROBLEM APPEARS UNSOLVABLE WITH GIVEN INPUT//) 
      GOTO 340 
!                                                                       
  220 WRITE (LOUT, 230) T, DTC 
  230 FORMAT(//35H KFLAG = -2 FROM INTEGRATOR AT T = ,E16.8,6H  DT =,   &
     &  E16.8/52H  THE REQUESTED ERROR IS SMALLER THAN CAN BE HANDLED//)
      GOTO 340 
!                                                                       
  240 WRITE (LOUT, 250) T 
  250 FORMAT(//37H INTEGRATION HALTED BY DRIVER AT T = ,E16.8/          &
     &   56H  EPS TOO SMALL TO BE ATTAINED FOR THE MACHINE PRECISION/)  
      KFLAG = - 2 
      GOTO 340 
!                                                                       
  260 WRITE (LOUT, 270) T 
  270 FORMAT(//35H KFLAG = -3 FROM INTEGRATOR AT T = ,E16.8/            &
     &   45H  CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED/)             
      GOTO 180 
!                                                                       
  280 WRITE (LOUT, 290) 
  290 FORMAT(//28H SINGULAR MATRIX ENCOUNTERED,                         &
     &         35H PROBABLY DUE TO STORAGE OVERWRITES//)                
      KFLAG = - 4 
      GOTO 340 
!                                                                       
  300 WRITE (LOUT, 310) T, TOUT, DTC 
  310 FORMAT(//45H INDEX = -1 ON INPUT WITH (T-TOUT)*DT .GE. 0./        &
     &   4H T =,E16.8,9H   TOUT =,E16.8,8H   DTC =,E16.8/               &
     &   44H INTERPOLATION WAS DONE AS ON NORMAL RETURN./               &
     &   41H DESIRED PARAMETER CHANGES WERE NOT MADE.)                  
      CALL INTERP (TOUT, WORK (IW10), NEQN, WORK (IW6) ) 
      INDEX = - 5 
      RETURN 
!                                                                       
  320 WRITE (LOUT, 330) IERID 
  330 FORMAT(//24H ILLEGAL INPUT...INDEX= ,I3//) 
      INDEX = IERID 
      RETURN 
!                                                                       
  335 WRITE (LOUT, 336) IWSTOR, IWSAVE, IISTOR, IISAVE 
  336 FORMAT(//21H INSUFFICIENT STORAGE/24H  WORK MUST BE OF LENGTH,    &
     & I10,5X,12HYOU PROVIDED,I10/24H IWORK MUST BE OF LENGTH,I10,5X,   &
     & 12HYOU PROVIDED,I10//)                                           
      INDEX = IERID 
      RETURN 
!                                                                       
  340 TOUT = T 
      DO 350 I = 1, NEQN 
        I1 = I - 1 
  350 WORK (IW6 + I1) = WORK (IW10 + I1) 
  360 INDEX = KFLAG 
      TOUTP = TOUT 
      DT = DTUSED 
      IF (KFLAG.NE.0) DT = DTC 
      RETURN 
      END SUBROUTINE PDECOL                         
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE VALUES (X, USOL, SCTCH, NDIM1, NDIM2, NPTS, NDERV,     &
      WORK)                                                             
      use FPwork
      SAVE 	 
!-----------------------------------------------------------------------
! SUBROUTINE VALUES COMPUTES THE SOLUTION U AND THE FIRST NDERV         
! DERIVATIVES OF U AT THE NPTS POINTS X AND AT TIME TOUT AND RETURNS    
! THEM IN THE ARRAY USOL.  THIS ROUTINE MUST BE USED TO OBTAIN          
! SOLUTION VALUES SINCE PDECOL DOES NOT RETURN ANY SOLUTION VALUES      
! TO THE USER.  SEE PDECOL.                                             
!                                                                       
! THE CALLING PARAMETERS ARE...                                         
!   X     =  AN ARBITRARY VECTOR OF SPATIAL POINTS OF LENGTH NPTS AT    
!            WHICH THE SOLUTION AND THE FIRST NDERV DERIVATIVE VALUES   
!            ARE TO BE CALCULATED.  IF X .LT. XLEFT ( X .GT. XRIGHT )   
!            THEN THE PIECEWISE POLYNOMIAL OVER THE LEFTMOST ( RIGHT-   
!            MOST ) INTERVAL IS EVALUATED TO CALCULATE THE SOLUTION     
!            VALUES AT THIS UNUSUAL VALUE OF  X.  SEE PDECOL.           
!                                                                       
!   USOL  =  AN ARRAY WHICH CONTAINS THE SOLUTION AND THE FIRST         
!            NDERV DERIVATIVES OF THE SOLUTION AT ALL THE POINTS IN     
!            THE INPUT VECTOR X.  IN PARTICULAR, USOL(J,I,K) CONTAINS   
!            THE VALUE OF THE (K-1)-ST DERIVATIVE OF THE J-TH PDE       
!            COMPONENT AT THE I-TH POINT OF THE X VECTOR FOR            
!            J = 1 TO NPDE, I = 1 TO NPTS, AND K = 1 TO NDERV+1.        
!                                                                       
!   SCTCH =  A USER SUPPLIED WORKING STORAGE ARRAY OF LENGTH AT LEAST   
!            KORD*(NDERV+1).  SEE BELOW AND PDECOL FOR DEFINITIONS OF   
!            THESE PARAMETERS.                                          
!                                                                       
!   NDIM1 =  THE FIRST DIMENSION OF THE OUTPUT ARRAY USOL IN THE CALLING
!            PROGRAM.  NDIM1 MUST BE .GE. NPDE.                         
!                                                                       
!   NDIM2 =  THE SECOND DIMENSION OF THE OUTPUT ARRAY USOL IN THE       
!            CALLING PROGRAM.  NDIM2 MUST BE .GE. NPTS.                 
!                                                                       
!   NPTS  =  THE NUMBER OF POINTS IN THE X VECTOR.                      
!                                                                       
!   NDERV =  THE NUMBER OF DERIVATIVE VALUES OF THE SOLUTION THAT ARE   
!            TO BE CALCULATED.  NDERV SHOULD BE LESS THAN KORD SINCE    
!            THE KORD-TH DERIVATIVE OF A POLYNOMIAL OF DEGREE KORD-1    
!            IS EQUAL TO ZERO.  SEE PDECOL.                             
!                                                                       
!   WORK  =  THE USERS WORKING STORAGE ARRAY WHICH IS USED IN THIS CASE 
!            TO PROVIDE THE CURRENT BASIS FUNCTION COEFFICIENTS AND THE 
!            PIECEWISE POLYNOMIAL BREAKPOINT SEQUENCE.                  
!                                                                       
! PACKAGE ROUTINES CALLED..  BSPLVD,INTERV                              
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                USERS MAIN PROGRAM                         
! FORTRAN FUNCTIONS USED..   NONE                                       
!                                                                       
!-----------------------------------------------------------------------
      DIMENSION USOL (NDIM1, NDIM2, NDERV), X (NPTS), SCTCH ( * ),      &
      WORK ( * )                                                        
      COMMON / SIZES / NINT, KORD, NCC, NPDE, NCPTS, NEQN, IQUAD 
!old      COMMON / ISTART / IW1, IW2, IW3, IW4, IW5, IW6, IDUM (12) 
      DATA ILEFT / 0 /, MFLAG / 0 / 
      NDERV1 = NDERV + 1 
      DO IPTS = 1, NPTS 
      CALL INTERV (WORK (IW2), NCPTS, X (IPTS), ILEFT, MFLAG) 
      CALL BSPLVD (WORK (IW2), KORD, X (IPTS), ILEFT, SCTCH, NDERV1) 
      IK = ILEFT - KORD 
      DO M = 1, NDERV1 
      I1 = (M - 1) * KORD 
      DO K = 1, NPDE 
      USOL (K, IPTS, M) = 0. 
      DO I = 1, KORD 
      I2 = (I + IK - 1) * NPDE+IW6 - 1 
      USOL (K, IPTS, M) = USOL (K, IPTS, M) + WORK (I2 + K) * SCTCH (I +&
      I1)                                                               
      enddo 
      enddo 
      enddo 
      enddo 
      RETURN 
      END SUBROUTINE VALUES                         
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE INITAL (K, A, RHS, X, XT, XC, PW, IPIV, ILEFT) 
           ! all variables                                              
      use FPuser  
      SAVE 	 
!-----------------------------------------------------------------------
! INITAL IS CALLED ONLY ONCE BY PDECOL TO PERFORM INITIALIZATION TASKS. 
! THESE TASKS INCLUDE - 1) DEFINING THE PIECEWISE POLYNOMIAL SPACE      
! BREAKPOINT SEQUENCE, 2) CALLING THE SUBROUTINE COLPNT TO DEFINE THE   
! REQUIRED COLLOCATION POINTS, 3) DEFING THE PIECEWISE POLYNOMIAL SPACE 
! BASIS FUNCTION VALUES (PLUS FIRST AND SECOND DERIVATIVE VALUES) AT    
! THE COLLOCATION POINTS, AND 4) DEFINING THE INITIAL BASIS FUNCTION    
! COEFFICIENTS WHICH DETERMINE THE PIECEWISE POLYNOMIAL WHICH           
! INTERPOLATES THE USER SUPPLIED (UINIT) INITIAL CONDITION FUNCTION(S)  
! AT THE COLLOCATION POINTS.                                            
!                                                                       
! K     = ORDER OF PIECEWISE POLYNOMIAL SPACE.                          
! A     = BASIS FUNCTION VALUES GENERATED BY INITAL.                    
! RHS   = TEMPORARY STORAGE USED TO RETURN INITIAL CONDITION COEFFICIENT
!         VALUES.                                                       
! X     = USER DEFINED PIECEWISE POLYNOMIAL BREAKPOINTS.                
! XT    = PIECEWISE POLYNOMIAL BREAKPOINT SEQUENCE GENERATED BY INITAL. 
! XC    = COLLOCATION POINTS GENERATED BY INITAL.                       
! PW    = STORAGE FOR BAND MATRIX USED TO GENERATE INITIAL              
!         COEFFICIENT VALUES.                                           
! IPIV  = PIVOT INFORMATION FOR LINEAR EQUATION SOLVER DECB-SOLB.       
! ILEFT = POINTERS TO BREAKPOINT SEQUENCE GENERATED BY INITAL.          
!                                                                       
! PACKAGE ROUTINES CALLED..  BSPLVD,COLPNT,DECB,INTERV,SOLB             
! USER ROUTINES CALLED..     UINIT                                      
! CALLED BY..                PDECOL                                     
! FORTRAN FUNCTIONS USED..   MAX0,MIN0                                  
!-----------------------------------------------------------------------
      DIMENSION A (K, 3, 1), RHS ( * ), X ( * ), XT ( * ), XC ( * ),    &
      PW ( * ), IPIV (8), ILEFT ( * )                                   
      COMMON / SIZES / NINT, KORD, NCC, NPDE, NCPTS, NEQN, IER 
      COMMON / GEAR9 / EPSJ, R0, ML, MU, IDUM (3), N0W 
      MFLAG = - 2 
      IER = 0 
!-----------------------------------------------------------------------
! SET UP THE PIECEWISE POLYNOMIAL SPACE BREAKPOINT SEQUENCE.            
!-----------------------------------------------------------------------
      KRPT = KORD-NCC 
      DO 10 I = 1, KORD 
        XT (NCPTS + I) = X (NINT + 1) 
   10 XT (I) = X (1) 
      DO 20 I = 2, NINT 
        I1 = (I - 2) * KRPT + KORD 
        DO 20 J = 1, KRPT 
   20 XT (I1 + J) = X (I) 
!-----------------------------------------------------------------------
! SET UP COLLOCATION POINTS ARRAY XC.                                   
!-----------------------------------------------------------------------
      CALL COLPNT (X, XC, XT) 
!-----------------------------------------------------------------------
! GENERATE THE ILEFT ARRAY.  STORE THE BASIS FUNCTION VALUES IN THE     
! ARRAY A.  THE ARRAY A IS DIMENSIONED A(KORD,3,NCPTS) AND A(K,J,I)     
! CONTAINS THE VALUE OF THE (J-1)-ST DERIVATIVE (J = 1,2,3) OF THE K-TH 
! NONZERO BASIS FUNCTION (K = 1, ... ,KORD) AT THE I-TH COLLOCATION     
! POINT (I = 1, ... ,NCPTS).  SET UP RHS FOR INTERPOLATING THE INITIAL  
! CONDITIONS AT THE COLLOCATION POINTS.  SET THE INTERPOLATION MATRIX   
! INTO THE BANDED MATRIX PW.                                            
!-----------------------------------------------------------------------
      DO 30 I = 1, N0W 
   30 PW (I) = 0. 
      DO 40 I = 1, NCPTS 
        CALL INTERV (XT, NCPTS, XC (I), ILEFT (I), MFLAG) 
        CALL BSPLVD (XT, KORD, XC (I), ILEFT (I), A (1, 1, I), 3) 
        I1 = NPDE * (I - 1) 
!old        CALL UINIT(XC(I),RHS(I1+1),NPDE)                            
        CALL UINIT (XC (I), RHS (I1 + 1) ) 
        ICOL = ILEFT (I) - I - 1 
        JL = MAX0 (1, I + 2 - NCPTS) 
        JU = MIN0 (KORD, KORD+I - 2) 
        DO 40 J = JL, JU 
          J1 = I1 + NEQN * (NPDE * (ICOL + J) - 1) 
          DO 40 JJ = 1, NPDE 
   40 PW (JJ + J1) = A (J, 1, I) 
!-----------------------------------------------------------------------
! LU DECOMPOSE THE MATRIX PW.                                           
!-----------------------------------------------------------------------
      CALL DECB (NEQN, NEQN, ML, MU, PW, IPIV, IER) 
      IF (IER.NE.0) RETURN 
!-----------------------------------------------------------------------
! SOLVE THE LINEAR SYSTEM   PW*Z = RHS.  THIS GIVES THE BASIS FUNCTION  
! COEFFICIENTS FOR THE INITIAL CONDITIONS.                              
!-----------------------------------------------------------------------
      CALL SOLB (NEQN, NEQN, ML, MU, PW, RHS, IPIV) 
      RETURN 
      END SUBROUTINE INITAL                         
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE COLPNT (X, XC, XT) 
        use FPglobal, ONLY: NOGAUS, MAXDER
        SAVE 	 
!-----------------------------------------------------------------------
! COLPNT IS CALLED ONLY ONCE BY INITAL TO DEFINE THE REQUIRED COLLOCA-  
! TION POINTS WHICH ARE TO BE USED WITH THE USER SELECTED PIECEWISE     
! POLYNOMIAL SPACE.  THE COLLOCATION POINTS ARE CHOSEN SUCH THAT THEY   
! ARE EITHER THE POINTS AT WHICH THE PIECEWISE POLYNOMIAL SPACE BASIS   
! FUNCTIONS ATTAIN THEIR UNIQUE MAXIMUM VALUES, OR, THE GAUSS-LEGENDRE  
! QUADRATURE POINTS WITHIN EACH PIECEWISE POLYNOMIAL SPACE SUBINTERVAL, 
! DEPENDING UPON THE SPACE BEING USED AND THE DESIRE OF THE USER.       
!                                                                       
! X  = USER DEFINED PIECEWISE POLYNOMIAL BREAKPOINTS.                   
! XC = COLLOCATION POINTS DEFINED BY COLPNT.                            
! XT = PIECEWISE POLYNOMIAL BREAKPOINT SEQUENCE.                        
!                                                                       
! PACKAGE ROUTINES CALLED..  BSPLVD,INTERV                              
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                INITAL                                     
! FORTRAN FUNCTIONS USED..   NONE                                       
!-----------------------------------------------------------------------
      DIMENSION RHO (40), X ( * ), XC ( * ), XT ( * ) 
      COMMON / SIZES / NINT, KORD, NCC, NPDE, NCPTS, NEQN, IQUAD 
!old      COMMON / OPTION / NOGAUS, MAXDER 
      DATA ILEFT / 0 / 
!-----------------------------------------------------------------------
! IF THE VARIABLE NOGAUS IN THE COMMON BLOCK /OPTION/ IS SET .EQ. 1,    
! THE USE OF THE GAUSS-LEGENDRE POINTS IS PROHIBITED FOR ALL CASES.     
! NOGAUS IS CURRENTLY SET .EQ. 0 BY A DATA STATEMENT IN THE BLOCK DATA. 
! THE USER MAY CHANGE THIS AS DESIRED.                                  
!-----------------------------------------------------------------------
      IF (NCC.NE.2.OR.NOGAUS.EQ.1) GOTO 200 
!-----------------------------------------------------------------------
! COMPUTE THE COLLOCATION POINTS TO BE AT THE GAUSS-LEGENDRE POINTS IN  
! EACH PIECEWISE POLYNOMIAL SPACE SUBINTERVAL.  THE ARRAY RHO IS SET TO 
! CONTAIN THE GAUSS-LEGENDRE POINTS FOR THE STANDARD INTERVAL (-1,1).   
!-----------------------------------------------------------------------
      IPTS = KORD-2 
      GOTO (10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,&
      150, 160, 170, 180), IPTS                                         
   10 RHO (1) = 0. 
      GOTO 190 
   20 RHO (2) = .577350269189626E-00 
      RHO (1) = - RHO (2) 
      GOTO 190 
   30 RHO (3) = .774596669241483E-00 
      RHO (1) = - RHO (3) 
      RHO (2) = 0. 
      GOTO 190 
   40 RHO (3) = .339981043584856E-00 
      RHO (2) = - RHO (3) 
      RHO (4) = .861136311594053E-00 
      RHO (1) = - RHO (4) 
      GOTO 190 
   50 RHO (4) = .538469310105683E-00 
      RHO (2) = - RHO (4) 
      RHO (5) = .906179845938664E-00 
      RHO (1) = - RHO (5) 
      RHO (3) = 0. 
      GOTO 190 
   60 RHO (4) = .238619186083197E-00 
      RHO (3) = - RHO (4) 
      RHO (5) = .661209386466265E-00 
      RHO (2) = - RHO (5) 
      RHO (6) = .932469514203152E-00 
      RHO (1) = - RHO (6) 
      GOTO 190 
   70 RHO (5) = .405845151377397E-00 
      RHO (3) = - RHO (5) 
      RHO (6) = .741531185599394E-00 
      RHO (2) = - RHO (6) 
      RHO (7) = .949107912342759E-00 
      RHO (1) = - RHO (7) 
      RHO (4) = 0. 
      GOTO 190 
   80 RHO (5) = .183434642495650E-00 
      RHO (4) = - RHO (5) 
      RHO (6) = .525532409916329E-00 
      RHO (3) = - RHO (6) 
      RHO (7) = .796666477413627E-00 
      RHO (2) = - RHO (7) 
      RHO (8) = .960289856497536E-00 
      RHO (1) = - RHO (8) 
      GOTO 190 
   90 RHO (5) = .0 
      RHO (6) = .324253423403809E-00 
      RHO (7) = .613371432700590E-00 
      RHO (8) = .836031107326636E-00 
      RHO (9) = .968160239507626E-00 
      DO 95 I = 1, 4 
   95 RHO (I) = - RHO (10 - I) 
      GOTO 190 
  100 RHO (6) = .148874338981631E-00 
      RHO (7) = .433395394129247E-00 
      RHO (8) = .679409568299024E-00 
      RHO (9) = .865063366688984E-00 
      RHO (10) = .973906528517172E-00 
      DO 105 I = 1, 5 
  105 RHO (I) = - RHO (11 - I) 
      GOTO 190 
  110 RHO (6) = .0 
      RHO (7) = .269543155952345E-00 
      RHO (8) = .519096129206812E-00 
      RHO (9) = .730152005574049E-00 
      RHO (10) = .887062599768095E-00 
      RHO (11) = .978228658146057E-00 
      DO 115 I = 1, 5 
  115 RHO (I) = - RHO (12 - I) 
      GOTO 190 
  120 RHO (7) = .125233408511469E-00 
      RHO (8) = .367831498998180E-00 
      RHO (9) = .587317954286617E-00 
      RHO (10) = .769902674194305E-00 
      RHO (11) = .904117256370475E-00 
      RHO (12) = .981560634246719E-00 
      DO 125 I = 1, 6 
  125 RHO (I) = - RHO (13 - I) 
      GOTO 190 
  130 RHO (7) = .0 
      RHO (8) = .230458315955135E-00 
      RHO (9) = .448492751036447E-00 
      RHO (10) = .642349339440340E-00 
      RHO (11) = .801578090733310E-00 
      RHO (12) = .917598399222978E-00 
      RHO (13) = .984183054718588E-00 
      DO 135 I = 1, 6 
  135 RHO (I) = - RHO (14 - I) 
      GOTO 190 
  140 RHO (8) = .108054948707344E-00 
      RHO (9) = .319112368927890E-00 
      RHO (10) = .515248636358154E-00 
      RHO (11) = .687292904811685E-00 
      RHO (12) = .827201315069765E-00 
      RHO (13) = .928434883663574E-00 
      RHO (14) = .986283808696812E-00 
      DO 145 I = 1, 7 
  145 RHO (I) = - RHO (15 - I) 
      GOTO 190 
  150 RHO (8) = .0 
      RHO (9) = .201194093997435E-00 
      RHO (10) = .394151347077563E-00 
      RHO (11) = .570972172608539E-00 
      RHO (12) = .724417731360170E-00 
      RHO (13) = .848206583410427E-00 
      RHO (14) = .937273392400706E-00 
      RHO (15) = .987992518020485E-00 
      DO 155 I = 1, 7 
  155 RHO (I) = - RHO (16 - I) 
      GOTO 190 
  160 RHO (9) = .950125098376374E-01 
      RHO (10) = .281603550779259E-00 
      RHO (11) = .458016777657227E-00 
      RHO (12) = .617876244402644E-00 
      RHO (13) = .755404408355003E-00 
      RHO (14) = .865631202387832E-00 
      RHO (15) = .944575023073233E-00 
      RHO (16) = .989400934991650E-00 
      DO 165 I = 1, 8 
  165 RHO (I) = - RHO (17 - I) 
      GOTO 190 
  170 RHO (9) = .0 
      RHO (10) = .178484181495848E-00 
      RHO (11) = .351231763453876E-00 
      RHO (12) = .512690537086477E-00 
      RHO (13) = .657671159216691E-00 
      RHO (14) = .781514003896801E-00 
      RHO (15) = .880239153726986E-00 
      RHO (16) = .950675521768768E-00 
      RHO (17) = .990575475314417E-00 
      DO 175 I = 1, 8 
  175 RHO (I) = - RHO (18 - I) 
      GOTO 190 
  180 RHO (10) = .847750130417353E-01 
      RHO (11) = .251886225691506E-00 
      RHO (12) = .411751161462843E-00 
      RHO (13) = .559770831073948E-00 
      RHO (14) = .691687043060353E-00 
      RHO (15) = .803704958972523E-00 
      RHO (16) = .892602466497556E-00 
      RHO (17) = .955823949571398E-00 
      RHO (18) = .991565168420931E-00 
      DO 185 I = 1, 9 
  185 RHO (I) = - RHO (19 - I) 
!-----------------------------------------------------------------------
! COMPUTE THE GAUSS-LEGENDRE COLLOCATION POINTS IN EACH SUBINTERVAL.    
!-----------------------------------------------------------------------
  190 DO 195 I = 1, NINT 
        FAC = (X (I + 1) - X (I) ) * .5 
        DO 195 J = 1, IPTS 
          KNOT = IPTS * (I - 1) + J + 1 
  195 XC (KNOT) = X (I) + FAC * (RHO (J) + 1.) 
      XC (1) = X (1) 
      XC (NCPTS) = X (NINT + 1) 
      RETURN 
!-----------------------------------------------------------------------
! COMPUTE THE COLLOCATION POINTS TO BE AT THE POINTS WHERE THE BASIS    
! FUNCTIONS ATTAIN THEIR MAXIMA.  A BISECTION METHOD IS USED TO FIND    
! THE POINTS TO MACHINE PRECISION.  THIS PROCESS COULD BE SPEEDED UP    
! BY USING A SECANT METHOD IF DESIRED.                                  
!-----------------------------------------------------------------------
  200 ITOP = NCPTS - 1 
      MFLAG = - 2 
      XC (1) = X (1) 
      XC (NCPTS) = X (NINT + 1) 
      DO 240 I = 2, ITOP 
        XOLD = 1.E+20 
        XL = XT (I) 
        XR = XT (I + KORD) 
  210   XNEW = .5 * (XL + XR) 
        IF (XOLD.EQ.XNEW) GOTO 240 
        CALL INTERV (XT, NCPTS, XNEW, ILEFT, MFLAG) 
        CALL BSPLVD (XT, KORD, XNEW, ILEFT, RHO, 2) 
        DO 220 J = 1, KORD 
          IF (I.EQ.J + ILEFT - KORD) GOTO 230 
  220   END DO 
  230   XVAL = RHO (KORD+J) 
        IF (XVAL.EQ.0.0) XR = XNEW 
        IF (XVAL.GT.0.0) XL = XNEW 
        IF (XVAL.LT.0.0) XR = XNEW 
        XOLD = XNEW 
        GOTO 210 
  240 XC (I) = XR 
      RETURN 
      END SUBROUTINE COLPNT                         
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE BSPLVD (XT, K, X, ILEFT, VNIKX, NDERIV) 
           ! all variables                                              
      SAVE 	 
!-----------------------------------------------------------------------
! THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE        
! EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.        
! SEE REFERENCE BELOW.                                                  
!                                                                       
! CALCULATES THE VALUE AND THE FIRST NDERIV-1 DERIVATIVES OF ALL        
! B-SPLINES WHICH DO NOT VANISH AT X.  THE ROUTINE FILLS THE TWO-       
! DIMENSIONAL ARRAY VNIKX(J,IDERIV), J=IDERIV, ... ,K WITH NONZERO      
! VALUES OF B-SPLINES OF ORDER K+1-IDERIV, IDERIV=NDERIV, ... ,1, BY    
! REPEATED CALLS TO BSPLVN.                                             
!                                                                       
! XT     = PIECEWISE POLYNOMIAL BREAKPOINT SEQUENCE.                    
! K      = ORDER OF THE PIECEWISE POLYNOMIAL SPACE.                     
! X      = POINT AT WHICH THE B-SPLINE IS TO BE EVALUATED.              
! ILEFT  = POINTER TO THE BREAKPOINT SEQUENCE.                          
! VNIKX  = TABLE OF B-SPLINE VALUES AND DERIVATIVES.                    
! NDERIV = DETERMINES NUMBER OF DERIVATIVES TO BE GENERATED.            
!                                                                       
! REFERENCE                                                             
!                                                                       
!    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.        
!      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.            
!                                                                       
! PACKAGE ROUTINES CALLED..  BSPLVN                                     
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                COLPNT,INITAL,VALUES                       
! FORTRAN FUNCTIONS USED..   FLOAT,MAX0                                 
!-----------------------------------------------------------------------
      DIMENSION XT (1), VNIKX (K, NDERIV) 
      DIMENSION A (20, 20) 
      KO = K + 1 - NDERIV 
      CALL BSPLVN (XT, KO, 1, X, ILEFT, VNIKX (NDERIV, NDERIV) ) 
      IF (NDERIV.LE.1) GOTO 120 
      IDERIV = NDERIV 
      DO 20 I = 2, NDERIV 
        IDERVM = IDERIV - 1 
        DO 10 J = IDERIV, K 
   10   VNIKX (J - 1, IDERVM) = VNIKX (J, IDERIV) 
        IDERIV = IDERVM 
        CALL BSPLVN (XT, 0, 2, X, ILEFT, VNIKX (IDERIV, IDERIV) ) 
   20 END DO 
      DO 40 I = 1, K 
        DO 30 J = 1, K 
   30   A (I, J) = 0. 
   40 A (I, I) = 1. 
      KMD = K 
      DO 110 M = 2, NDERIV 
        KMD = KMD-1 
        FKMD = FLOAT (KMD) 
        I = ILEFT 
        J = K 
   50   JM1 = J - 1 
        IPKMD = I + KMD 
        DIFF = XT (IPKMD) - XT (I) 
        IF (JM1.EQ.0) GOTO 80 
        IF (DIFF.EQ.0.) GOTO 70 
        DO 60 L = 1, J 
   60   A (L, J) = (A (L, J) - A (L, J - 1) ) / DIFF * FKMD 
   70   J = JM1 
        I = I - 1 
        GOTO 50 
   80   IF (DIFF.EQ.0.) GOTO 90 
        A (1, 1) = A (1, 1) / DIFF * FKMD 
   90   DO 110 I = 1, K 
          V = 0. 
          JLOW = MAX0 (I, M) 
          DO 100 J = JLOW, K 
  100     V = A (I, J) * VNIKX (J, M) + V 
  110 VNIKX (I, M) = V 
  120 RETURN 
      END SUBROUTINE BSPLVD                         
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE BSPLVN (XT, JHIGH, INDEX, X, ILEFT, VNIKX) 
           ! all variables                                              
      SAVE 	 
!-----------------------------------------------------------------------
! THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE        
! EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.        
! SEE REFERENCE BELOW.                                                  
!                                                                       
! CALCULATES THE VALUE OF ALL POSSIBLY NONZERO B-SPLINES AT THE         
! POINT X OF ORDER MAX(JHIGH,(J+1)(INDEX-1)) FOR THE BREAKPOINT SEQ-    
! UENCE XT.  ASSUMING THAT XT(ILEFT) .LE. X .LE. XT(ILEFT+1), THE ROUT- 
! INE RETURNS THE B-SPLINE VALUES IN THE ONE DIMENSIONAL ARRAY VNIKX.   
!                                                                       
! FOR DEFINITIONS OF CALLING ARGUMENTS SEE ABOVE AND BSPLVD.            
!                                                                       
! REFERENCE                                                             
!                                                                       
!    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.        
!      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.            
!                                                                       
! PACKAGE ROUTINES CALLED..  NONE                                       
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                BSPLVD                                     
! FORTRAN FUNCTIONS USED..   NONE                                       
!-----------------------------------------------------------------------
      DIMENSION XT (1), VNIKX (1) 
      DIMENSION DELTAM (20), DELTAP (20) 
      DATA J / 1 /, DELTAM / 20 * 0.E-00 /, DELTAP / 20 * 0.E-00 / 
      GOTO (10, 20), INDEX 
   10 J = 1 
      VNIKX (1) = 1. 
      IF (J.GE.JHIGH) GOTO 40 
   20 IPJ = ILEFT + J 
      DELTAP (J) = XT (IPJ) - X 
      IMJP1 = ILEFT - J + 1 
      DELTAM (J) = X - XT (IMJP1) 
      VMPREV = 0. 
      JP1 = J + 1 
      DO 30 L = 1, J 
        JP1ML = JP1 - L 
        VM = VNIKX (L) / (DELTAP (L) + DELTAM (JP1ML) ) 
        VNIKX (L) = VM * DELTAP (L) + VMPREV 
   30 VMPREV = VM * DELTAM (JP1ML) 
      VNIKX (JP1) = VMPREV 
      J = JP1 
      IF (J.LT.JHIGH) GOTO 20 
   40 RETURN 
      END SUBROUTINE BSPLVN                         
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE INTERV (XT, LXT, X, ILEFT, MFLAG) 
           ! all variables                                              
      SAVE 	 
!-----------------------------------------------------------------------
! THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE        
! EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.        
! SEE REFERENCE BELOW.                                                  
!                                                                       
! COMPUTES LARGEST ILEFT IN (1,LXT) SUCH THAT XT(ILEFT) .LE. X.  THE    
! PROGRAM STARTS THE SEARCH FOR ILEFT WITH THE VALUE OF ILEFT THAT WAS  
! RETURNED AT THE PREVIOUS CALL (AND WAS SAVED IN THE LOCAL VARIABLE    
! ILO) TO MINIMIZE THE WORK IN THE COMMON CASE THAT THE VALUE OF X ON   
! THIS CALL IS CLOSE TO THE VALUE OF X ON THE PREVIOUS CALL.  SHOULD    
! THIS ASSUMPTION NOT BE VALID, THEN THE PROGRAM LOCATES ILO AND IHI    
! SUCH THAT XT(ILO) .LE. X .LT. XT(IHI) AND, ONCE THEY ARE FOUND USES   
! BISECTION TO FIND THE CORRECT VALUE FOR ILEFT. MFLAG IS AN ERROR FLAG.
!                                                                       
! FOR DEFINITIONS OF CALLING ARGUMENTS SEE ABOVE AND BSPLVD.            
!                                                                       
! REFERENCE                                                             
!                                                                       
!    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.        
!      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.            
!                                                                       
! PACKAGE ROUTINES CALLED..  NONE                                       
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                COLPNT,INITAL,VALUES                       
! FORTRAN FUNCTIONS USED..   NONE                                       
!-----------------------------------------------------------------------
      DIMENSION XT (LXT) 
      IF (MFLAG.EQ. - 2) ILO = 1 
      IHI = ILO + 1 
      IF (IHI.LT.LXT) GOTO 20 
      IF (X.GE.XT (LXT) ) GOTO 110 
      IF (LXT.LE.1) GOTO 90 
      ILO = LXT - 1 
      GOTO 21 
   20 IF (X.GE.XT (IHI) ) GOTO 40 
   21 IF (X.GE.XT (ILO) ) GOTO 100 
!-----------------------------------------------------------------------
! NOW X .LT. XT(IHI).  FIND LOWER BOUND.                                
!-----------------------------------------------------------------------
   30 ISTEP = 1 
   31 IHI = ILO 
      ILO = IHI - ISTEP 
      IF (ILO.LE.1) GOTO 35 
      IF (X.GE.XT (ILO) ) GOTO 50 
      ISTEP = ISTEP * 2 
      GOTO 31 
   35 ILO = 1 
      IF (X.LT.XT (1) ) GOTO 90 
      GOTO 50 
!-----------------------------------------------------------------------
! NOW X .GE. XT(ILO).  FIND UPPER BOUND.                                
!-----------------------------------------------------------------------
   40 ISTEP = 1 
   41 ILO = IHI 
      IHI = ILO + ISTEP 
      IF (IHI.GE.LXT) GOTO 45 
      IF (X.LT.XT (IHI) ) GOTO 50 
      ISTEP = ISTEP * 2 
      GOTO 41 
   45 IF (X.GE.XT (LXT) ) GOTO 110 
      IHI = LXT 
!-----------------------------------------------------------------------
! NOW XT(ILO) .LE. X .LT. XT(IHI).  NARROW THE INTERVAL.                
!-----------------------------------------------------------------------
   50 MIDDLE = (ILO + IHI) / 2 
      IF (MIDDLE.EQ.ILO) GOTO 100 
!-----------------------------------------------------------------------
! NOTE..  IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1.          
!-----------------------------------------------------------------------
      IF (X.LT.XT (MIDDLE) ) GOTO 53 
      ILO = MIDDLE 
      GOTO 50 
   53 IHI = MIDDLE 
      GOTO 50 
!-----------------------------------------------------------------------
! SET OUTPUT AND RETURN.                                                
!-----------------------------------------------------------------------
   90 MFLAG = - 1 
      ILEFT = 1 
      RETURN 
  100 MFLAG = 0 
      ILEFT = ILO 
      RETURN 
  110 MFLAG = 1 
      ILEFT = LXT 
      RETURN 
      END SUBROUTINE INTERV                         
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE STIFIB (N0, Y, YMAX, ERROR, SAVE1, SAVE2, SAVE3, PW,   &
      IPIV, WORK, IWORK)                                                
        use FPglobal, ONLY: NOGAUS, MAXDER
        use FPwork
        SAVE 	 
!-----------------------------------------------------------------------
! STIFIB PERFORMS ONE STEP OF THE INTEGRATION OF AN INITIAL VALUE       
! PROBLEM FOR A SYSTEM OF ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM,  
!   A(Y,T)*(DY/DT) = G(Y,T),   WHERE   Y = (Y(1),Y(2), ... ,Y(N)).      
! STIFIB IS FOR USE WHEN THE MATRICES A AND DG/DY HAVE BANDED OR NEARLY 
! BANDED FORM.  THE DEPENDENCE OF A(Y,T) ON Y IS ASSUMED TO BE WEAK.    
!                                                                       
! REFERENCE                                                             
!                                                                       
!   HINDMARSH, A.C., PRELIMINARY DOCUMENTATION OF GEARIB.. SOLUTION     
!     OF IMPLICIT SYSTEMS OF ORDINARY DIFFERENTIAL EQUATIONS WITH       
!     BANDED JACOBIANS, LAWRENCE LIVERMORE LAB, UCID-30130, FEBRUARY    
!     1976.                                                             
!                                                                       
! COMMUNICATION WITH STIFIB IS DONE WITH THE FOLLOWING VARIABLES..      
!                                                                       
!   Y       AN N0 BY LMAX ARRAY CONTAINING THE DEPENDENT VARIABLES      
!             AND THEIR SCALED DERIVATIVES.  LMAX IS 13 FOR THE ADAMS   
!             METHODS AND 6 FOR THE GEAR METHODS.  LMAX - 1 = MAXDER    
!             IS THE MAXIMUM ORDER AVAILABLE.  SEE SUBROUTINE CSTCOL.   
!             Y(I,J+1) CONTAINS THE J-TH DERIVATIVE OF Y(I), SCALED BY  
!             H**J/FACTORIAL(J)  (J = 0,1,...,NQ).                      
!   N0      A CONSTANT INTEGER .GE. N, USED FOR DIMENSIONING PURPOSES.  
!   T       THE INDEPENDENT VARIABLE. T IS UPDATED ON EACH STEP TAKEN.  
!   H       THE STEP SIZE TO BE ATTEMPTED ON THE NEXT STEP.             
!             H IS ALTERED BY THE ERROR CONTROL ALGORITHM DURING THE    
!             PROBLEM.  H CAN BE EITHER POSITIVE OR NEGATIVE, BUT ITS   
!             SIGN MUST REMAIN CONSTANT THROUGHOUT THE PROBLEM.         
!   HMIN,   THE MINIMUM AND MAXIMUM ABSOLUTE VALUE OF THE STEP SIZE     
!    HMAX     TO BE USED FOR THE STEP.  THESE MAY BE CHANGED AT ANY     
!             TIME, BUT WILL NOT TAKE EFFECT UNTIL THE NEXT H CHANGE.   
!   EPS     THE RELATIVE ERROR BOUND.  SEE DESCRIPTION IN PDECOL.       
!   UROUND  THE UNIT ROUNDOFF OF THE MACHINE.                           
!   N       THE NUMBER OF FIRST-ORDER DIFFERENTIAL EQUATIONS.           
!   MF      THE METHOD FLAG.  SEE DESCRIPTION IN PDECOL.                
!   KFLAG   A COMPLETION CODE WITH THE FOLLOWING MEANINGS..             
!                     0  THE STEP WAS SUCCESFUL.                        
!                    -1  THE REQUESTED ERROR COULD NOT BE ACHIEVED      
!                          WITH ABS(H) = HMIN.                          
!                    -2  THE REQUESTED ERROR IS SMALLER THAN CAN        
!                          BE HANDLED FOR THIS PROBLEM.                 
!                    -3  CORRECTOR CONVERGENCE COULD NOT BE             
!                          ACHIEVED FOR ABS(H) = HMIN.                  
!                    -4  SINGULAR A-MATRIX ENCOUNTERED.                 
!             ON A RETURN WITH KFLAG NEGATIVE, THE VALUES OF T AND      
!             THE Y ARRAY ARE AS OF THE BEGINNING OF THE LAST           
!             STEP, AND H IS THE LAST STEP SIZE ATTEMPTED.              
!   JSTART  AN INTEGER USED ON INPUT AND OUTPUT.                        
!             ON INPUT, IT HAS THE FOLLOWING VALUES AND MEANINGS..      
!                     0  PERFORM THE FIRST STEP.                        
!                 .GT.0  TAKE A NEW STEP CONTINUING FROM THE LAST.      
!                 .LT.0  TAKE THE NEXT STEP WITH A NEW VALUE OF         
!                          H, EPS, N, AND/OR MF.                        
!             ON EXIT, JSTART IS NQ, THE CURRENT ORDER OF THE METHOD.   
!   YMAX    AN ARRAY OF N ELEMENTS WITH WHICH THE ESTIMATED LOCAL       
!             ERRORS IN Y ARE COMPARED.                                 
!   ERROR   AN ARRAY OF N ELEMENTS.  ERROR(I)/TQ(2) IS THE ESTIMATED    
!             ONE-STEP ERROR IN Y(I).                                   
!   SAVE1,SAVE2,SAVE3   THREE WORKING STORAGE ARRAYS, EACH OF LENGTH N. 
!   PW      A BLOCK OF LOCATIONS USED FOR THE CHORD ITERATION           
!             MATRIX.  SEE DESCRIPTION IN PDECOL.                       
!   IPIV    AN INTEGER ARRAY OF LENGTH N FOR PIVOT INFORMATION.         
!   ML,MU   THE LOWER AND UPPER HALF BANDWIDTHS, RESPECTIVELY, OF       
!             THE CHORD ITERATION MATRIX.  SEE DESCRIPTION IN PDECOL.   
!   WORK,IWORK   WORKING ARRAYS WHICH ARE USED TO PASS APPROPRIATE      
!                  ARRAYS TO OTHER SUBROUTINES.                         
!                                                                       
! PACKAGE ROUTINES CALLED..  CSTCOL,DIFFUN,PSETIB,RES,SOLB              
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                PDECOL                                     
! FORTRAN FUNCTIONS USED..   ABS,AMAX1,AMIN1,FLOAT                      
!-----------------------------------------------------------------------
      DIMENSION Y (N0, * ), YMAX (N0), ERROR (N0), SAVE1 (N0), SAVE2 (  &
      N0), SAVE3 (N0), PW (1), IPIV ( * ), WORK ( * ), IWORK ( * )      
      COMMON / SIZES / NINT, KORD, NCC, NPDE, NCPTS, NEQN, IQUAD 
!old      COMMON / ISTART / IW1, IW2, IW3, IW4, IW5, IW6, IW7, IW8, IW9,    &
!old      IW10, IW11, IW12, IW13, IW14, IW15, IW16, IW17, IW18              
      COMMON / GEAR1 / T, H, HMIN, HMAX, EPS, UROUND, N, MF, KFLAG,     &
      JSTART                                                            
      COMMON / GEAR9 / EPSJ, R0, ML, MU, MW, NM1, N0ML, N0W 
      COMMON / GEAR0 / HUSED, NQUSED, NSTEP, NFE, NJE 
!old      COMMON / OPTION / NOGAUS, MAXDER 
      DIMENSION EL (13), TQ (4) 
      DATA EL (2) / 1. /, OLDL0 / 1. /, TQ (1) / 0. /, IER / 0 / 
      KFLAG = 0 
      TOLD = T 
      IF (JSTART.GT.0) GOTO 200 
      IF (JSTART.NE.0) GOTO 120 
!-----------------------------------------------------------------------
! ON THE FIRST CALL, THE ORDER IS SET TO 1 AND THE INITIAL YDOT IS      
! CALCULATED.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE INCREASED    
! IN A SINGLE STEP.  IT IS INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL   
! INITIAL H, BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE            
! OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT 2     
! FOR THE NEXT INCREASE.                                                
!-----------------------------------------------------------------------
      NQ = 1 
      IER = 0 
      CALL DIFFUN (N, T, Y, SAVE1, IER, PW, IPIV, WORK, IWORK) 
      IF (IER.NE.0) GOTO 685 
      DO 110 I = 1, N 
  110 Y (I, 2) = H * SAVE1 (I) 
      METH = MF / 10 
      MITER = MF - 10 * METH 
      L = 2 
      IDOUB = 3 
      RMAX = 1.E+04 
      RC = 0. 
      CRATE = 1. 
      EPSOLD = EPS 
      HOLD = H 
      MFOLD = MF 
      NOLD = N 
      NSTEP = 0 
      NSTEPJ = 0 
      NFE = 0 
      NJE = 1 
      IRET = 3 
      GOTO 130 
!-----------------------------------------------------------------------
! IF THE CALLER HAS CHANGED METH, CSTCOL IS CALLED TO SET               
! THE COEFFICIENTS OF THE METHOD.  IF THE CALLER HAS CHANGED            
! N, EPS, OR METH, THE CONSTANTS E, EDN, EUP, AND BND MUST BE RESET.    
! E IS A COMPARISON FOR ERRORS OF THE CURRENT ORDER NQ. EUP IS          
! TO TEST FOR INCREASING THE ORDER, EDN FOR DECREASING THE ORDER.       
! BND IS USED TO TEST FOR CONVERGENCE OF THE CORRECTOR ITERATES.        
! IF THE CALLER HAS CHANGED H, Y MUST BE RESCALED.                      
! IF H OR METH HAS BEEN CHANGED, IDOUB IS RESET TO L + 1 TO PREVENT     
! FURTHER CHANGES IN H FOR THAT MANY STEPS.                             
!-----------------------------------------------------------------------
  120 IF (MF.EQ.MFOLD) GOTO 150 
      MEO = METH 
      MIO = MITER 
      METH = MF / 10 
      MITER = MF - 10 * METH 
      MFOLD = MF 
      IF (MITER.NE.MIO) IWEVAL = MITER 
      IF (METH.EQ.MEO) GOTO 150 
      IDOUB = L + 1 
      IRET = 1 
  130 CALL CSTCOL (METH, NQ, EL, TQ) 
      LMAX = MAXDER + 1 
      RC = RC * EL (1) / OLDL0 
      OLDL0 = EL (1) 
  140 FN = FLOAT (N) 
      EDN = FN * (TQ (1) * EPS) **2 
      E = FN * (TQ (2) * EPS) **2 
      EUP = FN * (TQ (3) * EPS) **2 
      BND = FN * (TQ (4) * EPS) **2 
      GOTO (160, 170, 200), IRET 
  150 IF ( (EPS.EQ.EPSOLD) .AND. (N.EQ.NOLD) ) GOTO 160 
      EPSOLD = EPS 
      NOLD = N 
      IRET = 1 
      GOTO 140 
  160 IF (H.EQ.HOLD) GOTO 200 
      RH = H / HOLD 
      H = HOLD 
      IREDO = 3 
      GOTO 175 
  170 RH = AMAX1 (RH, HMIN / ABS (H) ) 
  175 RH = AMIN1 (RH, HMAX / ABS (H), RMAX) 
      R1 = 1. 
      DO 180 J = 2, L 
        R1 = R1 * RH 
        DO 180 I = 1, N 
  180 Y (I, J) = Y (I, J) * R1 
      H = H * RH 
      RC = RC * RH 
      IDOUB = L + 1 
      IF (IREDO.EQ.0) GOTO 690 
!-----------------------------------------------------------------------
! THIS SECTION COMPUTES THE PREDICTED VALUES BY EFFECTIVELY             
! MULTIPLYING THE Y ARRAY BY THE PASCAL TRIANGLE MATRIX.                
! RC IS THE RATIO OF NEW TO OLD VALUES OF THE COEFFICIENT  H*EL(1).     
! WHEN RC DIFFERS FROM 1 BY MORE THAN 30 PERCENT, OR THE CALLER HAS     
! CHANGED MITER, IWEVAL IS SET TO MITER TO FORCE PW TO BE UPDATED.      
! IN ANY CASE, PW IS UPDATED AT LEAST EVERY 40-TH STEP.                 
! PW IS THE CHORD ITERATION MATRIX A - H*EL(1)*(DG/DY).                 
!-----------------------------------------------------------------------
  200 IF (ABS (RC - 1.) .GT.0.3) IWEVAL = MITER 
      IF (NSTEP.GE.NSTEPJ + 40) IWEVAL = MITER 
      T = T + H 
      DO 210 J1 = 1, NQ 
        DO 210 J2 = J1, NQ 
          J = (NQ + J1) - J2 
          DO 210 I = 1, N 
  210 Y (I, J) = Y (I, J) + Y (I, J + 1) 
!-----------------------------------------------------------------------
! UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS        
! MADE ON THE R.M.S. NORM OF EACH CORRECTION, USING BND, WHICH          
! IS DEPENDENT ON EPS.  THE SUM OF THE CORRECTIONS IS ACCUMULATED       
! IN THE VECTOR ERROR(I).  THE Y ARRAY IS NOT ALTERED IN THE CORRECTOR  
! LOOP.  THE UPDATED Y VECTOR IS STORED TEMPORARILY IN SAVE1.           
! THE UPDATED H*YDOT IS STORED IN SAVE2.                                
!-----------------------------------------------------------------------
  220 DO 230 I = 1, N 
        SAVE2 (I) = Y (I, 2) 
  230 ERROR (I) = 0. 
      M = 0 
      CALL RES (T, H, Y, SAVE2, SAVE3, NPDE, NCPTS, WORK (IW1), IWORK,  &
      WORK, WORK (IW14), WORK (IW15), WORK (IW16), WORK (IW3), WORK (   &
      IW9) )                                                            
      NFE = NFE+1 
      IF (IWEVAL.LE.0) GOTO 350 
!-----------------------------------------------------------------------
! IF INDICATED, THE MATRIX PW IS REEVALUATED BEFORE STARTING THE        
! CORRECTOR ITERATION.  IWEVAL IS SET TO 0 AS AN INDICATOR              
! THAT THIS HAS BEEN DONE.  PW IS COMPUTED AND PROCESSED IN PSETIB.     
!-----------------------------------------------------------------------
      IWEVAL = 0 
      RC = 1. 
      NJE = NJE+1 
      NSTEPJ = NSTEP 
      CON = - H * EL (1) 
!old      CALL PSETIB (Y, PW, N0, CON, MITER, IER, WORK(IW1), IWORK,    
!old     *     WORK(IW3),WORK(IW9),SAVE2,IPIV,YMAX,WORK(IW11),WORK(IW12)
!old     *     WORK(IW13),WORK(IW16),WORK(IW14),WORK(IW15),WORK,NPDE)   
      CALL PSETIB (Y, PW, N0, CON, MITER, IER, WORK (IW1), IWORK, WORK (&
      IW3), WORK (IW9), SAVE2, IPIV, YMAX, WORK (IW11), WORK (IW12),    &
      WORK (IW13), WORK (IW16), WORK (IW14), WORK (IW15), WORK)         
      IF (IER.NE.0) GOTO 420 
!-----------------------------------------------------------------------
! COMPUTE THE CORRECTOR ERROR, R SUB M, AND SOLVE THE LINEAR SYSTEM     
! WITH THAT AS RIGHT-HAND SIDE AND PW AS COEFFICIENT MATRIX,            
! USING THE LU DECOMPOSITION OF PW.                                     
!-----------------------------------------------------------------------
  350 CALL SOLB (N0, N, ML, MU, PW, SAVE3, IPIV) 
  370 D = 0. 
      DO 380 I = 1, N 
        ERROR (I) = ERROR (I) + SAVE3 (I) 
        D = D+ (SAVE3 (I) / YMAX (I) ) **2 
        SAVE1 (I) = Y (I, 1) + EL (1) * ERROR (I) 
  380 SAVE2 (I) = Y (I, 2) + ERROR (I) 
!-----------------------------------------------------------------------
! TEST FOR CONVERGENCE.  IF M.GT.0, AN ESTIMATE OF THE CONVERGENCE      
! RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.       
!-----------------------------------------------------------------------
  400 IF (M.NE.0) CRATE = AMAX1 (.9 * CRATE, D / D1) 
      IF ( (D * AMIN1 (1.E0, 2. * CRATE) ) .LE.BND) GOTO 450 
      D1 = D 
      M = M + 1 
      IF (M.EQ.3) GOTO 410 
      CALL RES (T, H, SAVE1, SAVE2, SAVE3, NPDE, NCPTS, WORK (IW1),     &
      IWORK, WORK, WORK (IW14), WORK (IW15), WORK (IW16), WORK (IW3),   &
      WORK (IW9) )                                                      
      GOTO 350 
!-----------------------------------------------------------------------
! THE CORRECTOR ITERATION FAILED TO CONVERGE IN 3 TRIES.                
! IF THE MATRIX PW IS NOT UP TO DATE, IT IS REEVALUATED FOR THE         
! NEXT TRY.  OTHERWISE THE Y ARRAY IS RETRACTED TO ITS VALUES           
! BEFORE PREDICTION, AND H IS REDUCED, IF POSSIBLE.  IF NOT, A          
! NO-CONVERGENCE EXIT IS TAKEN.                                         
!-----------------------------------------------------------------------
  410 NFE = NFE+2 
      IF (IWEVAL.EQ. - 1) GOTO 440 
  420 T = TOLD 
      RMAX = 2. 
      DO 430 J1 = 1, NQ 
        DO 430 J2 = J1, NQ 
          J = (NQ + J1) - J2 
          DO 430 I = 1, N 
  430 Y (I, J) = Y (I, J) - Y (I, J + 1) 
      IF (ABS (H) .LE.HMIN * 1.00001) GOTO 680 
      RH = .25 
      IREDO = 1 
      GOTO 170 
  440 IWEVAL = MITER 
      GOTO 220 
!-----------------------------------------------------------------------
! THE CORRECTOR HAS CONVERGED.  IWEVAL IS SET TO -1 TO SIGNAL           
! THAT PW MAY NEED UPDATING ON SUBSEQUENT STEPS.  THE ERROR TEST        
! IS MADE AND CONTROL PASSES TO STATEMENT 500 IF IT FAILS.              
!-----------------------------------------------------------------------
  450 IWEVAL = - 1 
      NFE = NFE+M 
      D = 0. 
      DO 460 I = 1, N 
  460 D = D+ (ERROR (I) / YMAX (I) ) **2 
      IF (D.GT.E) GOTO 500 
!-----------------------------------------------------------------------
! AFTER A SUCCESSFUL STEP, UPDATE THE Y ARRAY.                          
! CONSIDER CHANGING H IF IDOUB = 1.  OTHERWISE DECREASE IDOUB BY 1.     
! IF IDOUB IS THEN 1 AND NQ .LT. MAXDER, THEN ERROR IS SAVED FOR        
! USE IN A POSSIBLE ORDER INCREASE ON THE NEXT STEP.                    
! IF A CHANGE IN H IS CONSIDERED, AN INCREASE OR DECREASE IN ORDER      
! BY ONE IS CONSIDERED ALSO.  A CHANGE IN H IS MADE ONLY IF IT IS BY A  
! FACTOR OF AT LEAST 1.1.  IF NOT, IDOUB IS SET TO 10 TO PREVENT        
! TESTING FOR THAT MANY STEPS.                                          
!-----------------------------------------------------------------------
      KFLAG = 0 
      IREDO = 0 
      NSTEP = NSTEP + 1 
      HUSED = H 
      NQUSED = NQ 
      DO 470 J = 1, L 
        DO 470 I = 1, N 
  470 Y (I, J) = Y (I, J) + EL (J) * ERROR (I) 
      IF (IDOUB.EQ.1) GOTO 520 
      IDOUB = IDOUB - 1 
      IF (IDOUB.GT.1) GOTO 700 
      IF (L.EQ.LMAX) GOTO 700 
      DO 490 I = 1, N 
  490 Y (I, LMAX) = ERROR (I) 
      GOTO 700 
!-----------------------------------------------------------------------
! THE ERROR TEST FAILED.  KFLAG KEEPS TRACK OF MULTIPLE FAILURES.       
! RESTORE T AND THE Y ARRAY TO THEIR PREVIOUS VALUES, AND PREPARE       
! TO TRY THE STEP AGAIN.  COMPUTE THE OPTIMUM STEP SIZE FOR THIS OR     
! ONE LOWER ORDER.                                                      
!-----------------------------------------------------------------------
  500 KFLAG = KFLAG - 1 
      T = TOLD 
      DO 510 J1 = 1, NQ 
        DO 510 J2 = J1, NQ 
          J = (NQ + J1) - J2 
          DO 510 I = 1, N 
  510 Y (I, J) = Y (I, J) - Y (I, J + 1) 
      RMAX = 2. 
      IF (ABS (H) .LE.HMIN * 1.00001) GOTO 660 
      IF (KFLAG.LE. - 3) GOTO 640 
      IREDO = 2 
      PR3 = 1.E+20 
      GOTO 540 
!-----------------------------------------------------------------------
! REGARDLESS OF THE SUCCESS OR FAILURE OF THE STEP, FACTORS             
! PR1, PR2, AND PR3 ARE COMPUTED, BY WHICH H COULD BE DIVIDED           
! AT ORDER NQ - 1, ORDER NQ, OR ORDER NQ + 1, RESPECTIVELY.             
! IN THE CASE OF FAILURE, PR3 = 1.E20 TO AVOID AN ORDER INCREASE.       
! THE SMALLEST OF THESE IS DETERMINED AND THE NEW ORDER CHOSEN          
! ACCORDINGLY.  IF THE ORDER IS TO BE INCREASED, WE COMPUTE ONE         
! ADDITIONAL SCALED DERIVATIVE.                                         
!-----------------------------------------------------------------------
  520 PR3 = 1.E+20 
      IF (L.EQ.LMAX) GOTO 540 
      D1 = 0. 
      DO 530 I = 1, N 
  530 D1 = D1 + ( (ERROR (I) - Y (I, LMAX) ) / YMAX (I) ) **2 
      ENQ3 = .5 / FLOAT (L + 1) 
      PR3 = ( (D1 / EUP) **ENQ3) * 1.4 + 1.4E-06 
  540 ENQ2 = .5 / FLOAT (L) 
      PR2 = ( (D / E) **ENQ2) * 1.2 + 1.2E-06 
      PR1 = 1.E+20 
      IF (NQ.EQ.1) GOTO 560 
      D = 0. 
      DO 550 I = 1, N 
  550 D = D+ (Y (I, L) / YMAX (I) ) **2 
      ENQ1 = .5 / FLOAT (NQ) 
      PR1 = ( (D / EDN) **ENQ1) * 1.3 + 1.3E-06 
  560 IF (PR2.LE.PR3) GOTO 570 
      IF (PR3.LT.PR1) GOTO 590 
      GOTO 580 
  570 IF (PR2.GT.PR1) GOTO 580 
      NEWQ = NQ 
      RH = 1. / PR2 
      GOTO 620 
  580 NEWQ = NQ - 1 
      RH = 1. / PR1 
      GOTO 620 
  590 NEWQ = L 
      RH = 1. / PR3 
      IF (RH.LT.1.1) GOTO 610 
      DO 600 I = 1, N 
  600 Y (I, NEWQ + 1) = ERROR (I) * EL (L) / FLOAT (L) 
      GOTO 630 
  610 IDOUB = 10 
      GOTO 700 
  620 IF ( (KFLAG.EQ.0) .AND. (RH.LT.1.1) ) GOTO 610 
!-----------------------------------------------------------------------
! IF THERE IS A CHANGE OF ORDER, RESET NQ, L, AND THE COEFFICIENTS.     
! IN ANY CASE H IS RESET ACCORDING TO RH AND THE Y ARRAY IS RESCALED.   
! THEN EXIT FROM 690 IF THE STEP WAS OK, OR REDO THE STEP OTHERWISE.    
!-----------------------------------------------------------------------
      IF (NEWQ.EQ.NQ) GOTO 170 
  630 NQ = NEWQ 
      L = NQ + 1 
      IRET = 2 
      GOTO 130 
!-----------------------------------------------------------------------
! CONTROL REACHES THIS SECTION IF 3 OR MORE FAILURES HAVE OCCURED.      
! IT IS ASSUMED THAT THE DERIVATIVES THAT HAVE ACCUMULATED IN THE       
! Y ARRAY HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST              
! DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO 1.  THEN            
! H IS REDUCED BY A FACTOR OF 10, AND THE STEP IS RETRIED.              
! AFTER A TOTAL OF 7 FAILURES, AN EXIT IS TAKEN WITH KFLAG = -2.        
!-----------------------------------------------------------------------
  640 IF (KFLAG.EQ. - 7) GOTO 670 
      RH = .1 
      RH = AMAX1 (HMIN / ABS (H), RH) 
      H = H * RH 
      IER = 0 
      CALL DIFFUN (N, T, Y, SAVE1, IER, PW, IPIV, WORK, IWORK) 
      IF (IER.NE.0) GOTO 685 
      NJE = NJE+1 
      DO 650 I = 1, N 
  650 Y (I, 2) = H * SAVE1 (I) 
      IWEVAL = MITER 
      IDOUB = 10 
      IF (NQ.EQ.1) GOTO 200 
      NQ = 1 
      L = 2 
      IRET = 3 
      GOTO 130 
!-----------------------------------------------------------------------
! ALL RETURNS ARE MADE THROUGH THIS SECTION.  H IS SAVED IN HOLD        
! TO ALLOW THE CALLER TO CHANGE H ON THE NEXT STEP.                     
!-----------------------------------------------------------------------
  660 KFLAG = - 1 
      GOTO 700 
  670 KFLAG = - 2 
      GOTO 700 
  680 KFLAG = - 3 
      GOTO 700 
  685 KFLAG = - 4 
      GOTO 700 
  690 RMAX = 10. 
  700 HOLD = H 
      JSTART = NQ 
      RETURN 
      END SUBROUTINE STIFIB                         
!///////////////////////////////////////////////////////////////////////
!old      SUBROUTINE GFUN (T,C,UDOT,NPDE,NCPTS,A,BC,DBDU,DBDUX,DZDT,    
!old     *                 XC,UVAL,ILEFT)                               
      SUBROUTINE GFUN (T, C, UDOT, A, BC, DBDU, DBDUX, DZDT, XC, UVAL,  &
      ILEFT)                                                            
!-----------------------------------------------------------------------
! CALLING ARGUMENTS ARE DEFINED BELOW AND IN PDECOL.                    
!                                                                       
! SUBROUTINE GFUN COMPUTES THE FUNCTION UDOT=G(C,T), THE RIGHT-         
! HAND SIDE OF THE SEMI-DISCRETE APPROXIMATION TO THE ORIGINAL          
! SYSTEM OF PARTIAL DIFFERENTIAL EQUATIONS AND UPDATES THE BOUNDARY     
! CONDITION INFORMATION.                                                
!                                                                       
! PACKAGE ROUTINES CALLED..  EVAL                                       
! USER ROUTINES CALLED..     BNDRY,F                                    
! CALLED BY..                DIFFUN,PSETIB,RES                          
! FORTRAN FUNCTIONS USED..   NONE                                       
!-----------------------------------------------------------------------
!      include"FPglobal.f90" 
      use FPglobal
      use FPgrid
      use FPuser
      SAVE 	 
!???      DIMENSION C (NPDE, NCPTS), UDOT (NPDE, NCPTS) 
      DIMENSION C (NPDE, *), UDOT (NPDE, *) 
!???      DIMENSION A(*), BC(NPDE,NPDE,4), XC(*), UVAL(NPDE,3),  ILEFT(*)
      DIMENSION A(*), BC(NPDE,NPDE,*), XC(*), UVAL(NPDE,*),  ILEFT(*)
      DIMENSION DZDT (NPDE), DBDU (NPDE, NPDE), DBDUX (NPDE, NPDE) 
!old      COMMON /SIZES/ NINT0,KORD0,IDUM(4),IQUAD0                     
!-----------------------------------------------------------------------
!old(                                                                   
!      REAL Xg (NCPTS), Ug (NPDE, NCPTS), UXg (NPDE, NCPTS), UXXg (NPDE, &
!      NCPTS)                                                            
!      REAL UTg (NPDE, NCPTS), Qg (NPDE, NCPTS) 
!      INTEGER iXg 
!      COMMON / UGRID / Xg, Ug, UXg, UXXg, UTg, Qg, iXg 
!old)                                                                   
!-----------------------------------------------------------------------
      DO 10 K = 1, 4 
        DO 10 J = 1, NPDE 
          DO 10 I = 1, NPDE 
            BC (I, J, K) = 0.0 
   10 CONTINUE 
!-----------------------------------------------------------------------
!new(                                                                   
      DO I = 1, NCPTS 
      Xg (I) = XC (I) 
      CALL EVAL (I, NPDE, C, UVAL, A, ILEFT) 
      DO j = 1, NPDE 
      Ug (j, I) = UVAL (j, 1) 
      UXg (j, I) = UVAL (j, 2) 
      UXXg (j, I) = UVAL (j, 3) 
      enddo 
      enddo 
!tmpdbg(                                                                
!      open (33,status='unknown')                                       
!      write (33,'("# T = ",1p,e9.3)') T                                
!      write (33,'(1p,4(1x,e12.5))')                                    
!     >   (Xg(i),Ug(1,i),UXg(1,i),UXXg(1,i),i=1,NCPTS)                  
!      close (33)                                                       
!tmpdbg)                                                                
                                                                        
      CALL FGRID (0, T) 
!new)                                                                   
!-----------------------------------------------------------------------
                                                                        
!-----------------------------------------------------------------------
! UPDATE THE LEFT BOUNDARY VALUES.  SAVE LEFT BOUNDARY CONDITION        
! INFORMATION IN THE FIRST 2*NPDE*NPDE LOCATIONS OF BC.                 
!                                                                       
! NOTE.. UVAL(K,1) = U(K), UVAL(K,2) = UX(K), AND UVAL(K,3) = UXX(K).   
!-----------------------------------------------------------------------
      CALL EVAL (1, NPDE, C, UVAL, A, ILEFT) 
!old      CALL BNDRY(T,XC(1),UVAL,UVAL(1,2),DBDU,DBDUX,DZDT,NPDE)       
!old      CALL F(T,XC(1),UVAL,UVAL(1,2),UVAL(1,3),UDOT,NPDE)            
      CALL BNDRY (T, XC (1), UVAL, UVAL (1, 2), DBDU, DBDUX, DZDT) 
      CALL F (1, UDOT (1, 1) ) 
      ILIM = KORD+2 
      DO 30 K = 1, NPDE 
        BC (K, K, 1) = 1. 
        IF (DBDU (K, K) .EQ.0.0.AND.DBDUX (K, K) .EQ.0.0) GOTO 30 
        UDOT (K, 1) = DZDT (K) 
        DO 20 J = 1, NPDE 
          BC (K, J, 2) = A (ILIM) * DBDUX (K, J) 
          BC (K, J, 1) = DBDU (K, J) - BC (K, J, 2) 
   20   END DO 
   30 END DO 
!-----------------------------------------------------------------------
! MAIN LOOP TO FORM RIGHT SIDE OF ODES AT THE COLLOCATION POINTS.       
!-----------------------------------------------------------------------
      ILIM = NCPTS - 1 
      DO 40 I = 2, ILIM 
        CALL EVAL (I, NPDE, C, UVAL, A, ILEFT) 
!old        CALL F(T,XC(I),UVAL,UVAL(1,2),UVAL(1,3),UDOT(1,I),NPDE)     
        CALL F (I, UDOT (1, I) ) 
   40 END DO 
!-----------------------------------------------------------------------
! UPDATE THE RIGHT BOUNDARY VALUES.  SAVE THE RIGHT BOUNDARY CONDITION  
! INFORMATION IN THE LAST 2*NPDE*NPDE LOCATIONS IN BC.                  
!-----------------------------------------------------------------------
      CALL EVAL (NCPTS, NPDE, C, UVAL, A, ILEFT) 
!old      CALL F(T,XC(NCPTS),UVAL,UVAL(1,2),UVAL(1,3),UDOT(1,NCPTS),NPDE
!old      CALL BNDRY(T,XC(NCPTS),UVAL,UVAL(1,2),DBDU,DBDUX,DZDT,NPDE)   
      CALL F (NCPTS, UDOT (1, NCPTS) ) 
      CALL BNDRY (T, XC (NCPTS), UVAL, UVAL (1, 2), DBDU, DBDUX, DZDT) 
      ILIM = NCPTS * 3 * KORD-KORD-1 
      DO 60 K = 1, NPDE 
        BC (K, K, 4) = 1. 
        IF (DBDU (K, K) .EQ.0.0.AND.DBDUX (K, K) .EQ.0.0) GOTO 60 
        UDOT (K, NCPTS) = DZDT (K) 
        DO 50 J = 1, NPDE 
          BC (K, J, 3) = A (ILIM) * DBDUX (K, J) 
          BC (K, J, 4) = DBDU (K, J) - BC (K, J, 3) 
   50   END DO 
   60 END DO 
      RETURN 
      END SUBROUTINE GFUN                           
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE EVAL (ICPT, NPDE, C, UVAL, A, ILEFT) 
           ! all variables                                              
      SAVE 	 
!-----------------------------------------------------------------------
! CALLING ARGUMENTS ARE DEFINED BELOW AND IN PDECOL.                    
!                                                                       
! SUBROUTINE EVAL EVALUATES U(K), UX(K), AND UXX(K), K=1 TO NPDE,       
! AT THE COLLOCATION POINT WITH INDEX ICPT USING THE VALUES OF          
! THE BASIS FUNCTION COEFFICIENTS IN C AND THE BASIS FUNCTION VALUES    
! STORED IN A.  THE RESULTS ARE STORED IN UVAL AS FOLLOWS..             
! UVAL(K,1) = U(K), UVAL(K,2) = UX(K), AND UVAL(K,3) = UXX(K).          
!                                                                       
! PACKAGE ROUTINES CALLED..  NONE                                       
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                GFUN,PDECOL,PSETIB                         
! FORTRAN FUNCTIONS USED..   NONE                                       
!-----------------------------------------------------------------------
      DIMENSION C (NPDE, * ), UVAL (NPDE, 3), A ( * ), ILEFT ( * ) 
      COMMON / SIZES / NINT, KORD, IDUM (5) 
      IK = ILEFT (ICPT) - KORD 
      IC = 3 * KORD * (ICPT - 1) 
      DO 10 M = 1, 3 
        ICC = IC + KORD * (M - 1) 
        DO 10 J = 1, NPDE 
          UVAL (J, M) = 0. 
          DO 10 I = 1, KORD 
            UVAL (J, M) = UVAL (J, M) + C (J, I + IK) * A (I + ICC) 
   10 CONTINUE 
      RETURN 
      END SUBROUTINE EVAL                           
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE DIFFUN (N, T, Y, YDOT, IER, PW, IPIV, WORK, IWORK) 
        use FPwork
        SAVE 	 
!-----------------------------------------------------------------------
! CALLING ARGUMENTS ARE DEFINED BELOW AND IN PDECOL.                    
!                                                                       
! THIS ROUTINE COMPUTES YDOT = A(Y,T)**-1 * G(Y,T) BY USE OF            
! THE ROUTINES GFUN, ADDA, DECB, AND SOLB.                              
!                                                                       
! PACKAGE ROUTINES CALLED..  ADDA,DECB,GFUN,SOLB                        
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                STIFIB                                     
! FORTRAN FUNCTIONS USED..   NONE                                       
!-----------------------------------------------------------------------
      DIMENSION Y (N), YDOT (N), PW (*), IPIV (*), WORK (*), IWORK (*) 
      COMMON / GEAR9 / EPSJ, R0, ML, MU, MW, NM1, N0ML, N0W 
      COMMON / SIZES / NINT, KORD, NCC, NPDE, NCPTS, NEQN, IQUAD 
!old      COMMON / ISTART / IW1, IW2, IW3, IDUM (5), IW9, IW10, IW11, IW12, &
!old      IW13, IW14, IW15, IW16, IW17, IW18                                
!old      CALL GFUN (T, Y, YDOT, NPDE, NCPTS, WORK(IW1), WORK, WORK(IW14
!old     *           WORK(IW15), WORK(IW16), WORK(IW3), WORK(IW9), IWORK
      CALL GFUN (T, Y, YDOT, WORK(IW1), WORK, WORK(IW14), WORK(IW15),&
      WORK(IW16), WORK(IW3), WORK(IW9), IWORK)                       
      DO I = 1, N0W 
      PW (I) = 0. 
      enddo 
      N0 = NM1 + 1 
      CALL ADDA (PW, N0, WORK (IW1), IWORK, WORK, NPDE) 
      CALL DECB (N0, N, ML, MU, PW, IPIV, IER) 
      IF (IER.NE.0) RETURN 
      CALL SOLB (N0, N, ML, MU, PW, YDOT, IPIV) 
      RETURN 
      END SUBROUTINE DIFFUN                         
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE ADDA (PW, N0, A, ILEFT, BC, NPDE) 
           ! all variables                                              
      SAVE 	 
!-----------------------------------------------------------------------
! CALLING ARGUMENTS ARE DEFINED BELOW AND IN PDECOL AND STIFIB.         
!                                                                       
! SUBROUTINE ADDA ADDS THE MATRIX  A  TO THE MATRIX STORED IN  PW  IN   
! BAND FORM.  PW IS STORED BY DIAGONALS WITH THE LOWERMOST DIAGONAL     
! STORED IN THE FIRST COLUMN OF THE ARRAY.                              
!                                                                       
! PACKAGE ROUTINES CALLED..  NONE                                       
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                DIFFUN,PSETIB                              
! FORTRAN FUNCTIONS USED..   NONE                                       
!-----------------------------------------------------------------------
      DIMENSION PW (N0, 1), A (1), ILEFT (1), BC (NPDE, NPDE, 4) 
      COMMON / SIZES / NINT, KORD, NCC, NPD, NCPTS, NEQN, IQUAD 
!-----------------------------------------------------------------------
! ADD THE BOUNDARY CONDITION PORTIONS OF THE A MATRIX TO PW ( THE FIRST 
! AND LAST BLOCK ROWS).                                                 
!-----------------------------------------------------------------------
      ICOL = (ILEFT (1) + IQUAD-1) * NPDE 
      DO I = 1, NPDE 
      IBOT = NEQN - NPDE+I 
      DO J = 1, NPDE 
      IND = ICOL + J - I 
      PW (I, IND) = PW (I, IND) + BC (I, J, 1) 
      PW (I, IND+NPDE) = PW (I, IND+NPDE) + BC (I, J, 2) 
      PW (IBOT, IND-NPDE) = PW (IBOT, IND-NPDE) + BC (I, J, 3) 
      PW (IBOT, IND) = PW (IBOT, IND) + BC (I, J, 4) 
      enddo 
      enddo 
!-----------------------------------------------------------------------
! UPDATE THE REMAINING ROWS OF PW BY ADDING THE APPROPRIATE VALUES      
! IN A TO PW.                                                           
!-----------------------------------------------------------------------
      IND = NCPTS - 1 
      DO I = 2, IND 
      I1 = (I - 1) * NPDE 
      I2 = (I - 1) * KORD * 3 
      ICOL = ILEFT (I) - I + IQUAD-1 
      DO J = 1, KORD 
      J1 = (ICOL + J) * NPDE 
      J2 = I2 + J 
      DO JJ = 1, NPDE 
      PW (I1 + JJ, J1) = PW (I1 + JJ, J1) + A (J2) 
      enddo 
      enddo 
      enddo 
      RETURN 
      END SUBROUTINE ADDA                           
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE RES (T, H, C, V, R, NPDE, NCPTS, A, ILEFT, BC, DBDU,   &
      DBDUX, DZDT, XC, UVAL)                                            
           ! all variables                                              
      SAVE 	 
!-----------------------------------------------------------------------
! CALLING ARGUMENTS ARE DEFINED BELOW AND IN PDECOL.                    
!                                                                       
! SUBROUTINE RES COMPUTES THE RESIDUAL VECTOR R = H*G(C,T) - A(C,T)*V   
! WHERE H IS THE CURRENT TIME STEP SIZE, G IS A VECTOR, A IS A          
! MATRIX, V IS A VECTOR, AND T IS THE CURRENT TIME.                     
!                                                                       
! PACKAGE ROUTINES CALLED..  GFUN                                       
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                STIFIB                                     
! FORTRAN FUNCTIONS USED..   NONE                                       
!-----------------------------------------------------------------------
      DIMENSION C (NPDE, NCPTS), R (NPDE, * ), V (NPDE, * ) 
      DIMENSION A ( * ), ILEFT ( * ), BC (NPDE, NPDE, 4), XC ( * ),     &
      UVAL ( * )                                                        
      DIMENSION DBDU (NPDE, NPDE), DBDUX (NPDE, NPDE), DZDT (NPDE) 
      COMMON / SIZES / NINT, KORD, NCC, IDUM (3), IQUAD 
!-----------------------------------------------------------------------
! FORM G(C,T) AND STORE IN R.                                           
!-----------------------------------------------------------------------
!old      CALL GFUN(T,C,R,NPDE,NCPTS,A,BC,DBDU,DBDUX,DZDT,XC,UVAL,ILEFT)
      CALL GFUN (T, C, R, A, BC, DBDU, DBDUX, DZDT, XC, UVAL, ILEFT) 
!-----------------------------------------------------------------------
! FORM THE FIRST AND LAST BLOCK ROWS OF THE RESIDUAL VECTOR             
! WHICH ARE DEPENDENT ON THE BOUNDARY CONDITIONS.                       
!-----------------------------------------------------------------------
      ILIM = NCPTS - 1 
      DO 20 I = 1, NPDE 
        SUM1 = 0.0 
        SUM2 = 0.0 
        DO 10 J = 1, NPDE 
          SUM1 = SUM1 + BC (I, J, 1) * V (J, 1) + BC (I, J, 2) * V (J,  &
          2)                                                            
          SUM2 = SUM2 + BC (I, J, 3) * V (J, ILIM) + BC (I, J, 4)       &
          * V (J, NCPTS)                                                
   10   END DO 
        R (I, 1) = H * R (I, 1) - SUM1 
        R (I, NCPTS) = H * R (I, NCPTS) - SUM2 
   20 END DO 
!-----------------------------------------------------------------------
! FORM THE REMAINING COMPONENTS OF THE RESIDUAL VECTOR.                 
!-----------------------------------------------------------------------
      DO 50 ICPTS = 2, ILIM 
        I2 = (ICPTS - 1) * KORD * 3 
        ICOL = ILEFT (ICPTS) - KORD 
        DO 40 JJ = 1, NPDE 
          SUM1 = 0. 
          DO 30 J = 1, KORD 
            SUM1 = SUM1 + A (I2 + J) * V (JJ, ICOL + J) 
   30     END DO 
          R (JJ, ICPTS) = H * R (JJ, ICPTS) - SUM1 
   40   END DO 
   50 END DO 
      RETURN 
      END SUBROUTINE RES                            
!///////////////////////////////////////////////////////////////////////
!old      SUBROUTINE PSETIB (C, PW, N0, CON, MITER, IER, A, ILEFT, XC, U
!old     *    SAVE2,IPIV,CMAX,DFDU,DFDUX,DFDUXX,DZDT,DBDU,DBDUX,BC,NPDE)
      SUBROUTINE PSETIB (C, PW, N0, CON, MITER, IER, A, ILEFT, XC, UVAL,&
      SAVE2, IPIV, CMAX, DFDU, DFDUX, DFDUXX, DZDT, DBDU, DBDUX, BC)    
!-----------------------------------------------------------------------
! CALLING ARGUMENTS ARE DEFINED BELOW AND IN PDECOL AND STIFIB.         
!                                                                       
! PSETIB IS CALLED BY STIFIB TO COMPUTE AND PROCESS THE MATRIX          
! PW = A - H*EL(1)*(DG/DC), WHERE A AND DG/DC ARE TREATED IN BAND       
! FORM.  DG/DC IS COMPUTED, EITHER WITH THE AID OF THE USER-SUPPLIED    
! ROUTINE DERIVF IF MITER = 1, OR BY FINITE DIFFERENCING WITH THE AID   
! OF THE PACKAGE-SUPPLIED ROUTINE DIFFF IF MITER = 2.  FINALLY,         
! PW IS SUBJECTED TO LU DECOMPOSITION IN PREPARATION FOR LATER          
! SOLUTION OF LINEAR SYSTEMS WITH PW AS COEFFICIENT MATRIX.             
! SEE SUBROUTINES DECB AND SOLB.                                        
!                                                                       
! IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION          
! WITH PSETIB USES THE FOLLOWING..                                      
!   EPSJ    = SQRT(UROUND), USED IN THE NUMERICAL JACOBIAN INCREMENTS.  
!   MW      = ML + MU + 1.                                              
!   NM1     = N0 - 1.                                                   
!   N0ML    = N0*ML.                                                    
!   N0W     = N0*MW.                                                    
!                                                                       
! PACKAGE ROUTINES CALLED..  ADDA,DECB,DIFFF,EVAL,GFUN                  
! USER ROUTINES CALLED..     BNDRY,DERIVF                               
! CALLED BY..                STIFIB                                     
! FORTRAN FUNCTIONS USED..   ABS,FLOAT,MAX0,MIN0,SQRT                   
!-----------------------------------------------------------------------
!      include"FPglobal.f90" 
      use FPglobal
      use FPgrid
      use FPuser
      SAVE 	 
      DIMENSION PW (N0, 1), C (1), CMAX (1) 
      DIMENSION A (1), ILEFT (1), BC (1), XC (1), UVAL (NPDE, 3),       &
      SAVE2 (1), IPIV (1)                                               
      DIMENSION DFDU (NPDE, NPDE), DFDUX (NPDE, NPDE), DFDUXX (NPDE,    &
      NPDE)                                                             
      DIMENSION DZDT (NPDE), DBDU (NPDE, NPDE), DBDUX (NPDE, NPDE) 
!old      COMMON /SIZES/ NINT,KORD,NCC,NPD,NCPTS,NEQN,IQUAD             
      COMMON / SIZES / idum1, idum2, idum3, NPD, idum4, NEQN, idum5 
      COMMON / GEAR1 / T, H, DUMMY (3), UROUND, N, IDUMMY (3) 
!old      COMMON /GEAR9/ EPSJ,R0,ML,MU,MW,NM1,N0ML,N0W                  
      COMMON / GEAR9 / EPSJ, R0, idum6, MU, MW, NM1, N0ML, N0W 
!old(                                                                   
!      REAL Xg (NCPTS), Ug (NPDE, NCPTS), UXg (NPDE, NCPTS), UXXg (NPDE, &
!      NCPTS)                                                            
!      REAL UTg (NPDE, NCPTS), Qg (NPDE, NCPTS) 
!      INTEGER iXg 
!      COMMON / UGRID / Xg, Ug, UXg, UXXg, UTg, Qg, iXg 
!old)                                                                   
!                                                                       
      DO I = 1, N0W 
      PW (I, 1) = 0. 
      enddo 
      IF (MITER.EQ.1) GOTO 25 
!old      CALL GFUN (T, C, SAVE2, NPDE, NCPTS,A,BC,DBDU,DBDUX,DZDT,XC,  
!old     *           UVAL,ILEFT)                                        
      CALL GFUN (T, C, SAVE2, A, BC, DBDU, DBDUX, DZDT, XC, UVAL, ILEFT) 
      D = 0. 
      DO I = 1, N 
      D = D+SAVE2 (I) **2 
      enddo 
      R0 = ABS (H) * SQRT (D / FLOAT (N0) ) * 1.E+03 * UROUND 
!-----------------------------------------------------------------------
! COMPUTE BLOCK ROWS OF JACOBIAN.                                       
!-----------------------------------------------------------------------
   25 CONTINUE 
!new(                                                                   
      DO I = 1, NCPTS 
      Xg (I) = XC (I) 
      CALL EVAL (I, NPDE, C, UVAL, A, ILEFT) 
      DO j = 1, NPDE 
      Ug (j, I) = UVAL (j, 1) 
      UXg (j, I) = UVAL (j, 2) 
      UXXg (j, I) = UVAL (j, 3) 
      enddo 
      enddo 
      CALL FGRID (0, T) 
!new)                                                                   
      DO 30 I = 1, NCPTS 
        I1 = (I - 1) * NPDE 
        I2 = (I - 1) * KORD * 3 
!old        CALL EVAL(I,NPDE,C,UVAL,A,ILEFT)                            
        IF (MITER.EQ.1) CALL DERIVF (T, XC (I), Ug (1, I), UXg (1, I),  &
        UXXg (1, I), DFDU, DFDUX, DFDUXX, NPDE)                         
!old     *      CALL DERIVF(T,XC(I),UVAL(1,1),UVAL(1,2),UVAL(1,3),      
!old     *                  DFDU,DFDUX,DFDUXX,NPDE)                     
!dbg DANGER!!!: test also for NPDE>1                                    
        IF (MITER.EQ.2) CALL DIFFF (T, XC (I), I, Ug (1, I), UXg (1, I),&
        UXXg (1, 3), DFDU, DFDUX, DFDUXX, CMAX, SAVE2)                  
!older     *       CALL DIFFF(T,XC(I),I,UVAL,UVAL(1,2),UVAL(1,3),       
!older     *                  DFDU,DFDUX,DFDUXX,NPDE,CMAX,SAVE2)        
!old     *       CALL DIFFF(T,XC(I),I,UVAL(1,1),UVAL(1,2),UVAL(1,3),    
!old     *                  DFDU,DFDUX,DFDUXX,CMAX,SAVE2)               
!dbg DANGER!!!: test also for NPDE>1                                    
        ICOL = ILEFT (I) - I + IQUAD-1 
        KLOW = MAX0 (1, I + 2 - NCPTS) 
        KUP = MIN0 (KORD, KORD+I - 2) 
        DO 30 KBLK = KLOW, KUP 
          J1 = (ICOL + KBLK) * NPDE 
          J2 = I2 + KBLK 
          J3 = J2 + KORD 
          J4 = J3 + KORD 
          DO 30 L = 1, NPDE 
            DO 30 K = 1, NPDE 
              PW (I1 + K, J1 - K + L) = DFDU (K, L) * A (J2) + DFDUX (K,&
              L) * A (J3) + DFDUXX (K, L) * A (J4)                      
   30 CONTINUE 
!-----------------------------------------------------------------------
! MODIFY THE LAST AND THE FIRST BLOCK ROWS FOR THE BOUNDARY CONDITIONS. 
! CURRENT INFORMATION FOR THE RIGHT BOUNDARY CONDITION IS ALREADY IN    
! THE ARRAYS DBDU, DBDUX AS A RESULT OF A PREVIOUS CALL TO GFUN.        
!-----------------------------------------------------------------------
      IROW = NEQN - NPDE 
      DO 50 K = 1, NPDE 
        IROW = IROW + 1 
        IF (DBDU (K, K) .EQ.0.0.AND.DBDUX (K, K) .EQ.0.0) GOTO 50 
        DO 40 J = 1, MW 
          PW (IROW, J) = 0.0 
   40   END DO 
   50 END DO 
      CALL EVAL (1, NPDE, C, UVAL, A, ILEFT) 
!old      CALL BNDRY(T,XC(1),UVAL,UVAL(1,2),DBDU,DBDUX,DZDT,NPDE)       
      CALL BNDRY (T, XC (1), UVAL, UVAL (1, 2), DBDU, DBDUX, DZDT) 
      DO 70 K = 1, NPDE 
        IF (DBDU (K, K) .EQ.0.0.AND.DBDUX (K, K) .EQ.0.0) GOTO 70 
        DO 60 J = 1, MW 
          PW (K, J) = 0.0 
   60   END DO 
   70 END DO 
      DO I = 1, N0W 
      PW (I, 1) = PW (I, 1) * CON 
      enddo 
!-----------------------------------------------------------------------
! ADD MATRIX A(C,T) TO PW.                                              
!-----------------------------------------------------------------------
      CALL ADDA (PW, N0, A, ILEFT, BC, NPDE) 
!-----------------------------------------------------------------------
! DO LU DECOMPOSITION ON PW.                                            
!-----------------------------------------------------------------------
      CALL DECB (N0, N, ML, MU, PW, IPIV, IER) 
      RETURN 
      END SUBROUTINE PSETIB                         
!///////////////////////////////////////////////////////////////////////
!old      SUBROUTINE DIFFF(T,X,IPT,U,UX,UXX,DFDU,DFDUX,DFDUXX,NPDE,CMAX,
!old     *                 SAVE2)                                       
      SUBROUTINE DIFFF (T, X, IPT, U, UX, UXX, DFDU, DFDUX, DFDUXX,     &
      CMAX, SAVE2)                                                      
!-----------------------------------------------------------------------
! CALLING ARGUMENTS ARE DEFINED BELOW AND IN PDECOL.                    
!                                                                       
! SUBROUTINE DIFFF IS USED IF MITER=2 TO PROVIDE FINITE DIFFERENCE      
! APPROXIMATIONS FOR THE PARTIAL DERIVATIVES OF THE K-TH USER DEFINED   
! FUNCTION IN THE F ROUTINE WITH RESPECT TO THE VARIABLES U, UX, AND    
! UXX.  THESE PARTIALS WITH RESPECT TO U, UX, AND UXX ARE COMPUTED,     
! STORED, AND RETURNED IN THE NPDE BY NPDE ARRAYS DFDU, DFDUX, AND      
! DFDUXX, RESPECTIVELY, AT COLLOCATION POINT NUMBER IPT.                
!                                                                       
! PACKAGE ROUTINES CALLED..  NONE                                       
! USER ROUTINES CALLED..     F                                          
! CALLED BY..                PSETIB                                     
! FORTRAN FUNCTIONS USED..   AMAX1                                      
!-----------------------------------------------------------------------
!      include "FPglobal.f90" 
      use FPglobal
      use FPgrid
      use FPuser
      SAVE 	 
      REAL U (NPDE), UX (NPDE), UXX (NPDE), DFDU (NPDE, NPDE), DFDUX (  &
      NPDE, NPDE), DFDUXX (NPDE, NPDE), CMAX (1), SAVE2 (1)             
      COMMON / GEAR9 / EPSJ, R0, ML0, MU, MW, NM1, N0ML, N0W 
!old(                                                                   
!      REAL Xg (NCPTS), Ug (NPDE, NCPTS), UXg (NPDE, NCPTS), UXXg (NPDE, &
!      NCPTS)                                                            
!      REAL UTg (NPDE, NCPTS), Qg (NPDE, NCPTS) 
!      INTEGER iXg 
!      COMMON / UGRID / Xg, Ug, UXg, UXXg, UTg, Qg, iXg 
!old)                                                                   
      ID = (IPT - 1) * NPDE 
      DO J = 1, NPDE 
      R = EPSJ * CMAX (J) 
      R = AMAX1 (R, R0) 
      RINV = 1. / R 
!                                                                       
!old(                                                                   
!        UJ = U(J)                                                      
!        U(J) = U(J) + R                                                
!old)                                                                   
!new(                                                                   
      UJ = Ug (J, IPT) 
      Ug (J, IPT) = Ug (J, IPT) + R 
      CALL FGRID (IPT, T) 
!new)                                                                   
!old         CALL F(T,X,U,UX,UXX,DFDU(1,J),NPDE)                        
      CALL F (IPT, DFDU (1, J) ) 
      DO I = 1, NPDE 
      DFDU (I, J) = (DFDU (I, J) - SAVE2 (I + ID) ) * RINV 
      enddo 
!old(                                                                   
!        U(J) = UJ                                                      
!        UJ = UX(J)                                                     
!        UX(J) = UX(J) + R                                              
!old)                                                                   
!new(                                                                   
      Ug (J, IPT) = UJ 
      UJ = UXg (J, IPT) 
      UXG (J, IPT) = UXG (J, IPT) + R 
      CALL FGRID (IPT, T) 
!new)                                                                   
!old         CALL F(T,X,U,UX,UXX,DFDUX(1,J),NPDE)                       
      CALL F (IPT, DFDUX (1, J) ) 
      DO I = 1, NPDE 
      DFDUX (I, J) = (DFDUX (I, J) - SAVE2 (I + ID) ) * RINV 
      enddo 
!old(                                                                   
!        UX(J) = UJ                                                     
!        UJ = UXX(J)                                                    
!        UXX(J) = UXX(J) + R                                            
!old)                                                                   
!new(                                                                   
      UXG (J, IPT) = UJ 
      UJ = UXXg (J, IPT) 
      UXXg (J, IPT) = UXXg (J, IPT) + R 
      CALL FGRID (IPT, T) 
!new)                                                                   
!old         CALL F(T,X,U,UX,UXX,DFDUXX(1,J),NPDE)                      
      CALL F (IPT, DFDUXX (1, J) ) 
      DO I = 1, NPDE 
      DFDUXX (I, J) = (DFDUXX (I, J) - SAVE2 (I + ID) ) * RINV 
      enddo 
!old(                                                                   
!       UXX(J) = UJ                                                     
!old)                                                                   
!new(                                                                   
      UXXg (J, IPT) = UJ 
!new)                                                                   
      enddo 
      RETURN 
      END SUBROUTINE DIFFF                          
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE INTERP (TOUT, Y, N0, Y0) 
           ! all variables                                              
      SAVE 	 
      DIMENSION Y0 (N0), Y (N0, 1) 
      COMMON / GEAR1 / T, H, DUMMY (4), N, IDUMMY (2), JSTART 
!-----------------------------------------------------------------------
!  SUBROUTINE INTERP COMPUTES INTERPOLATED VALUES OF THE DEPENDENT      
!  VARIABLE Y AND STORES THEM IN Y0.  THE INTERPOLATION IS TO THE       
!  POINT T = TOUT, AND USES THE NORDSIECK HISTORY ARRAY Y, AS FOLLOWS.. 
!                              NQ                                       
!                   Y0(I)  =  SUM  Y(I,J+1)*S**J ,                      
!                             J=0                                       
!  WHERE S = -(T-TOUT)/H.                                               
!-----------------------------------------------------------------------
      DO 10 I = 1, N 
   10 Y0 (I) = Y (I, 1) 
      L = JSTART + 1 
      S = (TOUT - T) / H 
      S1 = 1. 
      DO 30 J = 2, L 
        S1 = S1 * S 
        DO 20 I = 1, N 
   20   Y0 (I) = Y0 (I) + S1 * Y (I, J) 
   30 END DO 
      RETURN 
!----------------------- END OF SUBROUTINE INTERP ----------------------
      END SUBROUTINE INTERP                         
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE CSTCOL (METH, NQ, EL, TQ) 
           ! all variables                                              
      SAVE 	 
!-----------------------------------------------------------------------
! CSTCOL IS CALLED BY THE INTEGRATOR AND SETS COEFFICIENTS USED THERE.  
! THE VECTOR EL, OF LENGTH NQ + 1, DETERMINES THE BASIC METHOD.         
! THE VECTOR TQ, OF LENGTH 4, IS INVOLVED IN ADJUSTING THE STEP SIZE    
! IN RELATION TO TRUNCATION ERROR.  ITS VALUES ARE GIVEN BY THE         
! PERTST ARRAY.                                                         
!                                                                       
! THE VECTORS EL AND TQ DEPEND ON METH AND NQ.                          
! THE MAXIMUM ORDER, MAXDER, OF THE METHODS AVAILABLE IS CURRENTLY      
! 12 FOR THE ADAMS METHODS AND 5 FOR THE BDF METHODS.  MAXDER DEFAULTS  
! TO 5 UNLESS THE USER SETS MAXDER TO SOME OTHER LEGITIMATE VALUE       
! THROUGH THE COMMON BLOCK /OPTION/.  SEE PDECOL FOR ADDITIONAL DETAILS.
! LMAX = MAXDER + 1 IS THE NUMBER OF COLUMNS IN THE Y ARRAY (SEE STIFIB 
! AND THE VARIABLE C, Y, OR WORK(IW10) IN PDECOL.                       
!                                                                       
! THE COEFFICIENTS IN PERTST NEED BE GIVEN TO ONLY ABOUT                
! ONE PERCENT ACCURACY.  THE ORDER IN WHICH THE GROUPS APPEAR BELOW     
! IS..  COEFFICIENTS FOR ORDER NQ - 1, COEFFICIENTS FOR ORDER NQ,       
! COEFFICIENTS FOR ORDER NQ + 1.  WITHIN EACH GROUP ARE THE             
! COEFFICIENTS FOR THE ADAMS METHODS, FOLLOWED BY THOSE FOR THE         
! BDF METHODS.                                                          
!                                                                       
! REFERENCE                                                             
!                                                                       
!   GEAR, C.W., NUMERICAL INITIAL VALUE PROBLEMS IN ORDINARY            
!     DIFFERENTIAL EQUATIONS, PRENTICE-HALL, ENGLEWOOD CLIFFS,          
!     N. J., 1971.                                                      
!                                                                       
! PACKAGE ROUTINES CALLED..  NONE                                       
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                STIFIB                                     
! FORTRAN FUNCTIONS USED..   FLOAT                                      
!-----------------------------------------------------------------------
      DIMENSION PERTST (12, 2, 3), EL (13), TQ (4) 
      DATA PERTST / 1., 1., 2., 1., .3158, .07407, .01391, .002182,     &
      .0002945, .00003492, .000003692, .0000003524, 1., 1., .5, .1667,  &
      .04167, 1., 1., 1., 1., 1., 1., 1., 2., 12., 24., 37.89, 53.33,   &
      70.08, 87.97, 106.9, 126.7, 147.4, 168.8, 191.0, 2.0, 4.5, 7.333, &
      10.42, 13.7, 1., 1., 1., 1., 1., 1., 1., 12.0, 24.0, 37.89, 53.33,&
      70.08, 87.97, 106.9, 126.7, 147.4, 168.8, 191.0, 1., 3.0, 6.0,    &
      9.167, 12.5, 1., 1., 1., 1., 1., 1., 1., 1. /                     
!                                                                       
      GOTO (1, 2), METH 
    1 GOTO (101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112),&
      NQ                                                                
    2 GOTO (201, 202, 203, 204, 205), NQ 
!-----------------------------------------------------------------------
! THE FOLLOWING COEFFICIENTS SHOULD BE DEFINED TO MACHINE ACCURACY.     
! FOR A GIVEN ORDER NQ, THEY CAN BE CALCULATED BY USE OF THE            
! GENERATING POLYNOMIAL L(T), WHOSE COEFFICIENTS ARE EL(I)..            
!      L(T) = EL(1) + EL(2)*T + ... + EL(NQ+1)*T**NQ.                   
! FOR THE IMPLICIT ADAMS METHODS, L(T) IS GIVEN BY                      
!      DL/DT = (T+1)*(T+2)* ... *(T+NQ-1)/K,    L(-1) = 0,              
! WHERE                 K = FACTORIAL(NQ-1).                            
! FOR THE BDF METHODS,                                                  
!      L(T) = (T+1)*(T+2)* ... *(T+NQ)/K,                               
! WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).               
!                                                                       
! THE ORDER IN WHICH THE GROUPS APPEAR BELOW IS..                       
! IMPLICIT ADAMS METHODS OF ORDERS 1 TO 12,                             
! BDF METHODS OF ORDERS 1 TO 5.                                         
!-----------------------------------------------------------------------
  101 EL (1) = 1.0E-00 
      GOTO 900 
  102 EL (1) = 0.5E-00 
      EL (3) = 0.5E-00 
      GOTO 900 
  103 EL (1) = 4.1666666666667E-01 
      EL (3) = 0.75E-00 
      EL (4) = 1.6666666666667E-01 
      GOTO 900 
  104 EL (1) = 0.375E-00 
      EL (3) = 9.1666666666667E-01 
      EL (4) = 3.3333333333333E-01 
      EL (5) = 4.1666666666667E-02 
      GOTO 900 
  105 EL (1) = 3.4861111111111E-01 
      EL (3) = 1.0416666666667E-00 
      EL (4) = 4.8611111111111E-01 
      EL (5) = 1.0416666666667E-01 
      EL (6) = 8.3333333333333E-03 
      GOTO 900 
  106 EL (1) = 3.2986111111111E-01 
      EL (3) = 1.1416666666667E-00 
      EL (4) = 0.625E-00 
      EL (5) = 1.7708333333333E-01 
      EL (6) = 0.025E-00 
      EL (7) = 1.3888888888889E-03 
      GOTO 900 
  107 EL (1) = 3.1559193121693E-01 
      EL (3) = 1.225E-00 
      EL (4) = 7.5185185185185E-01 
      EL (5) = 2.5520833333333E-01 
      EL (6) = 4.8611111111111E-02 
      EL (7) = 4.8611111111111E-03 
      EL (8) = 1.9841269841270E-04 
      GOTO 900 
  108 EL (1) = 3.0422453703704E-01 
      EL (3) = 1.2964285714286E-00 
      EL (4) = 8.6851851851852E-01 
      EL (5) = 3.3576388888889E-01 
      EL (6) = 7.7777777777778E-02 
      EL (7) = 1.0648148148148E-02 
      EL (8) = 7.9365079365079E-04 
      EL (9) = 2.4801587301587E-05 
      GOTO 900 
  109 EL (1) = 2.9486800044092E-01 
      EL (3) = 1.3589285714286E-00 
      EL (4) = 9.7655423280423E-01 
      EL (5) = 0.4171875E-00 
      EL (6) = 1.1135416666667E-01 
      EL (7) = 0.01875E-00 
      EL (8) = 1.9345238095238E-03 
      EL (9) = 1.1160714285714E-04 
      EL (10) = 2.7557319223986E-06 
      GOTO 900 
  110 EL (1) = 2.8697544642857E-01 
      EL (3) = 1.4144841269841E-00 
      EL (4) = 1.0772156084656E-00 
      EL (5) = 4.9856701940035E-01 
      EL (6) = 0.1484375E-00 
      EL (7) = 2.9060570987654E-02 
      EL (8) = 3.7202380952381E-03 
      EL (9) = 2.9968584656085E-04 
      EL (10) = 1.3778659611993E-05 
      EL (11) = 2.7557319223986E-07 
      GOTO 900 
  111 EL (1) = 2.8018959644394E-01 
      EL (3) = 1.4644841269841E-00 
      EL (4) = 1.1715145502646E-00 
      EL (5) = 5.7935819003527E-01 
      EL (6) = 1.8832286155203E-01 
      EL (7) = 4.1430362654321E-02 
      EL (8) = 6.2111441798942E-03 
      EL (9) = 6.2520667989418E-04 
      EL (10) = 4.0417401528513E-05 
      EL (11) = 1.5156525573192E-06 
      EL (12) = 2.5052108385442E-08 
      GOTO 900 
  112 EL (1) = 2.7426554003160E-01 
      EL (3) = 1.5099386724387E-00 
      EL (4) = 1.2602711640212E-00 
      EL (5) = 6.5923418209877E-01 
      EL (6) = 2.3045800264550E-01 
      EL (7) = 5.5697246105232E-02 
      EL (8) = 9.4394841269841E-03 
      EL (9) = 1.1192749669312E-03 
      EL (10) = 9.0939153439153E-05 
      EL (11) = 4.8225308641975E-06 
      EL (12) = 1.5031265031265E-07 
      EL (13) = 2.0876756987868E-09 
      GOTO 900 
  201 EL (1) = 1.0E-00 
      GOTO 900 
  202 EL (1) = 6.6666666666667E-01 
      EL (3) = 3.3333333333333E-01 
      GOTO 900 
  203 EL (1) = 5.4545454545455E-01 
      EL (3) = EL (1) 
      EL (4) = 9.0909090909091E-02 
      GOTO 900 
  204 EL (1) = 0.48E-00 
      EL (3) = 0.7E-00 
      EL (4) = 0.2E-00 
      EL (5) = 0.02E-00 
      GOTO 900 
  205 EL (1) = 4.3795620437956E-01 
      EL (3) = 8.2116788321168E-01 
      EL (4) = 3.1021897810219E-01 
      EL (5) = 5.4744525547445E-02 
      EL (6) = 3.6496350364964E-03 
!                                                                       
  900 DO 910 K = 1, 3 
  910 TQ (K) = PERTST (NQ, METH, K) 
      TQ (4) = .5E-00 * TQ (2) / FLOAT (NQ + 2) 
      RETURN 
      END SUBROUTINE CSTCOL                         
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE DECB (NDIM, N, ML, MU, B, IPIV, IER) 
           ! all variables                                              
      SAVE 	 
!-----------------------------------------------------------------------
! SUBROUTINES DECB AND SOLB FORM A TWO SUBROUTINE PACKAGE FOR THE       
! DIRECT SOLUTION OF A SYSTEM OF LINEAR EQUATIONS IN WHICH THE          
! COEFFICIENT MATRIX IS REAL AND BANDED.                                
!                                                                       
!    LU DECOMPOSITION OF BAND MATRIX A..  L*U = P*A , WHERE P IS A      
!       PERMUTATION MATRIX, L IS A UNIT LOWER TRIANGULAR MATRIX,        
!       AND U IS AN UPPER TRIANGULAR MATRIX.                            
!    N     =  ORDER OF MATRIX.                                          
!    B     =  N BY (2*ML+MU+1) ARRAY CONTAINING THE MATRIX A ON INPUT   
!             AND ITS FACTORED FORM ON OUTPUT.                          
!             ON INPUT, B(I,K) (1.LE.I.LE.N) CONTAINS THE K-TH          
!             DIAGONAL OF A, OR A(I,J) IS STORED IN B(I,J-I+ML+1).      
!             ON OUTPUT, B CONTAINS THE L AND U FACTORS, WITH           
!             U IN COLUMNS 1 TO ML+MU+1, AND L IN COLUMNS               
!             ML+MU+2 TO 2*ML+MU+1.                                     
!    ML,MU =  WIDTHS OF THE LOWER AND UPPER PARTS OF THE BAND, NOT      
!             COUNTING THE MAIN DIAGONAL. TOTAL BANDWIDTH IS ML+MU+1.   
!    NDIM  =  THE FIRST DIMENSION (COLUMN LENGTH) OF THE ARRAY B.       
!             NDIM MUST BE .GE. N.                                      
!    IPIV  =  ARRAY OF LENGTH N CONTAINING PIVOT INFORMATION.           
!    IER   =  ERROR INDICATOR..                                         
!          =  0  IF NO ERROR,                                           
!          =  K  IF THE K-TH PIVOT CHOSEN WAS ZERO (A IS SINGULAR).     
!    THE INPUT ARGUMENTS ARE  NDIM, N, ML, MU, B.                       
!    THE OUTPUT ARGUMENTS ARE  B, IPIV, IER.                            
!                                                                       
!-----------------------------------------------------------------------
      DIMENSION B (NDIM, 1), IPIV (N) 
      IER = 0 
      IF (N.EQ.1) GOTO 92 
      LL = ML + MU + 1 
      N1 = N - 1 
      IF (ML.EQ.0) GOTO 32 
      DO 30 I = 1, ML 
        II = MU + I 
        K = ML + 1 - I 
        DO 10 J = 1, II 
   10   B (I, J) = B (I, J + K) 
        K = II + 1 
        DO 20 J = K, LL 
   20   B (I, J) = 0. 
   30 END DO 
   32 LR = ML 
      DO 90 NR = 1, N1 
        NP = NR + 1 
        IF (LR.NE.N) LR = LR + 1 
        MX = NR 
        XM = ABS (B (NR, 1) ) 
        IF (ML.EQ.0) GOTO 42 
        DO 40 I = NP, LR 
          IF (ABS (B (I, 1) ) .LE.XM) GOTO 40 
          MX = I 
          XM = ABS (B (I, 1) ) 
   40   END DO 
   42   IPIV (NR) = MX 
        IF (MX.EQ.NR) GOTO 60 
        DO 50 I = 1, LL 
          XX = B (NR, I) 
          B (NR, I) = B (MX, I) 
   50   B (MX, I) = XX 
   60   XM = B (NR, 1) 
        IF (XM.EQ.0.) GOTO 100 
        B (NR, 1) = 1. / XM 
        IF (ML.EQ.0) GOTO 90 
        XM = - B (NR, 1) 
        KK = MIN0 (N - NR, LL - 1) 
        DO 80 I = NP, LR 
          J = LL + I - NR 
          XX = B (I, 1) * XM 
          B (NR, J) = XX 
          DO 70 II = 1, KK 
   70     B (I, II) = B (I, II + 1) + XX * B (NR, II + 1) 
   80   B (I, LL) = 0. 
   90 END DO 
   92 NR = N 
      IF (B (N, 1) .EQ.0.) GOTO 100 
      B (N, 1) = 1. / B (N, 1) 
      RETURN 
  100 IER = NR 
      RETURN 
      END SUBROUTINE DECB                           
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE SOLB (NDIM, N, ML, MU, B, Y, IPIV) 
           ! all variables                                              
      SAVE 	 
!-----------------------------------------------------------------------
! SUBROUTINES DECB AND SOLB FORM A TWO SUBROUTINE PACKAGE FOR THE       
! DIRECT SOLUTION OF A SYSTEM OF LINEAR EQUATIONS IN WHICH THE          
! COEFFICIENT MATRIX IS REAL AND BANDED.                                
!                                                                       
!    SOLUTION OF  A*X = C  GIVEN LU DECOMPOSITION OF A FROM DECB.       
!    Y  =  RIGHT-HAND VECTOR C, OF LENGTH N, ON INPUT,                  
!       =  SOLUTION VECTOR X ON OUTPUT.                                 
!    ALL THE ARGUMENTS ARE INPUT ARGUMENTS.                             
!    THE OUTPUT ARGUMENT IS  Y.                                         
!                                                                       
! PACKAGE ROUTINES CALLED..  NONE                                       
! USER ROUTINES CALLED..     NONE                                       
! CALLED BY..                DIFFUN,INITAL,STIFIB                       
! FORTRAN FUNCTIONS USED..   MIN0                                       
!-----------------------------------------------------------------------
      DIMENSION B (NDIM, 1), Y (N), IPIV (N) 
      IF (N.EQ.1) GOTO 60 
      N1 = N - 1 
      LL = ML + MU + 1 
      IF (ML.EQ.0) GOTO 32 
      DO 30 NR = 1, N1 
        IF (IPIV (NR) .EQ.NR) GOTO 10 
        J = IPIV (NR) 
        XX = Y (NR) 
        Y (NR) = Y (J) 
        Y (J) = XX 
   10   KK = MIN0 (N - NR, ML) 
        DO 20 I = 1, KK 
   20   Y (NR + I) = Y (NR + I) + Y (NR) * B (NR, LL + I) 
   30 END DO 
   32 LL = LL - 1 
      Y (N) = Y (N) * B (N, 1) 
      KK = 0 
      DO 50 NB = 1, N1 
        NR = N - NB 
        IF (KK.NE.LL) KK = KK + 1 
        DP = 0. 
        IF (LL.EQ.0) GOTO 50 
        DO 40 I = 1, KK 
   40   DP = DP + B (NR, I + 1) * Y (NR + I) 
   50 Y (NR) = (Y (NR) - DP) * B (NR, 1) 
      RETURN 
   60 Y (1) = Y (1) * B (1, 1) 
      RETURN 
      END SUBROUTINE SOLB                           
!///////////////////////////////////////////////////////////////////////
!     f90 (self probing) realization of r1mach                          
      REAL function r1mach (i) 
      IMPLICIT none 
      INTEGER i 
      REAL b, x 
!                                                                       
!   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.           
!   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.           
!   R1MACH(3) = B**(-T), the smallest relative spacing.                 
!   R1MACH(4) = B**(1-T), the largest relative spacing.                 
!   R1MACH(5) = LOG10(B)                                                
!                                                                       
      x = 1.0 
      b = radix (x) 
      selectcase (i) 
      case (1) 
                                         ! the smallest positive magnitu
      r1mach = b** (minexponent (x) - 1) 
      case (2) 
                                         ! the largest magnitude.       
      r1mach = huge (x) 
      case (3) 
                                         ! the smallest relative spacing
      r1mach = b** ( - digits (x) ) 
      case (4) 
                                         ! the largest relative spacing.
      r1mach = b** (1 - digits (x) ) 
      case (5) 
      r1mach = log10 (b) 
      case default 
      WRITE ( * , '("ERROR in r1mach: i = ",i2," is out of bounds")') 
      STOP 
      endselect 
      RETURN 
      END FUNCTION r1mach                           
end module FPpdecol
