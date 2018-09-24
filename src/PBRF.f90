
SUBROUTINE KM(N,Y1,Y2,D,NT,TIM,LT,LC,FUN)
   IMPLICIT NONE
   INTEGER,INTENT(IN)::N,NT,D(N),LT,LC
   REAL(8),INTENT(IN)::Y1(N),Y2(N),TIM(NT)
   REAL(8),INTENT(OUT)::FUN(NT)
   INTEGER::R(N),I
   REAL(8)::LAM(N)
   IF (LT==0) THEN
      DO I=1,N,1
         R(I)=COUNT(Y2>=Y2(I))
	  END DO
   ELSE
      DO I=1,N,1
         R(I)=COUNT(Y2>=Y2(I) .AND. Y2(I)>=Y1)
	  END DO
   END IF
   LAM=1.0-REAL(D,8)/REAL(R,8)

   IF (LC==0) THEN
      DO I=1,NT,1
         FUN(I)=PRODUCT(LAM,Y2<=TIM(I))
      END DO
   ELSE
      DO I=1,NT,1
         FUN(I)=PRODUCT(LAM,Y2<TIM(I))
      END DO
   END IF
   RETURN
END SUBROUTINE KM


SUBROUTINE YPBF(N,Y1,Y2,D1,D2,NT,TIM,PBRF,V3,V4)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'ypbf_' ::YPBF
   IMPLICIT NONE
   INTEGER,INTENT(IN)::N,NT,D1(N),D2(N)
   REAL(8),INTENT(IN)::Y1(N),Y2(N),TIM(NT)
   REAL(8),INTENT(OUT)::PBRF(NT,4),V3(NT,3),V4(NT,3)

   REAL(8)::ST(NT),STY(N),G1(N),G2(N),STYG1(N),STYG2(N)
   REAL(8),ALLOCATABLE::Y11(:),Y12(:)
   INTEGER,ALLOCATABLE::DD(:)
   INTEGER::I,J,K,U,N1

   REAL(8)::VPBRF3(N,NT,3),VPBRF4(N,NT,3),RR(N),W1(N),W2(N)


   STYG1=0.0;STYG2=0.0
   CALL KM(N,Y1,Y1,1-D1,N,Y1,0,1,G1)
   CALL KM(N,Y2,Y2,1-D2,N,Y1,0,1,G2)

   N1=SUM(D1)
   ALLOCATE(Y11(N1),Y12(N1),DD(N1))
   Y11=PACK(Y1,D1==1)
   Y12=PACK(Y2,D1==1)
   DD=PACK(D2,D1==1)
   CALL KM(N1,Y11,Y12,DD,N,Y1,1,0,STY)
   CALL KM(N1,Y11,Y12,DD,NT,TIM,1,0,ST)
   DO I=1,N,1
      RR(I)=REAL(COUNT(Y11<=Y2(I) .AND. Y2(I)<=Y12),8)/REAL(N,8)
	  IF (RR(I)>0.0) THEN
	     RR(I)=1.0/RR(I)
	  END IF
	  IF (STY(I)>0.0 .AND. G1(I)>0.0) THEN
	     STYG1(I)=1.0/STY(I)/G1(I)
	  END IF
	  IF (STY(I)>0.0 .AND. G2(I)>0.0) THEN
	     STYG2(I)=1.0/STY(I)/G2(I)
	  END IF
	  W1(I)=REAL(COUNT(Y1>=Y1(I)),8)/REAL(N,8)
	  W2(I)=REAL(COUNT(Y2>=Y2(I)),8)/REAL(N,8)
   END DO
   DEALLOCATE(Y11,Y12,DD)

   PBRF=0.0;V3=1.0;V4=1.0
   VPBRF3=0.0;VPBRF4=0.0
   DO U=1,NT,1
      PBRF(U,3)=ST(U)/REAL(N,8)*SUM(STYG1,D1==1 .AND. Y1<=TIM(U))
      PBRF(U,4)=ST(U)/REAL(N,8)*SUM(STYG2,D1==1 .AND. Y1<=TIM(U))
	  VPBRF3(:,U,1)=-PBRF(U,3)
      VPBRF4(:,U,1)=-PBRF(U,4)
      DO I=1,N,1
	     IF (D1(I)==1 .AND. Y1(I)<=TIM(U)) THEN
	        VPBRF3(I,U,1)=VPBRF3(I,U,1)+ST(U)*STYG1(I)
            VPBRF4(I,U,1)=VPBRF4(I,U,1)+ST(U)*STYG2(I)
		 END IF

		 IF (D1(I)==1 .AND. D2(I)==1 .AND. Y2(I)<=TIM(U)) THEN
		    VPBRF3(I,U,2)=VPBRF3(I,U,2)-ST(U)*RR(I)*SUM(STYG1,D1==1 .AND. Y1<Y2(I))/REAL(N,8)
			VPBRF4(I,U,2)=VPBRF4(I,U,2)-ST(U)*RR(I)*SUM(STYG2,D1==1 .AND. Y1<Y2(I))/REAL(N,8)
         END IF

		 IF (D1(I)==1) THEN
            DO J=1,N,1
			   IF (Y1(I)<=Y2(J) .AND. Y2(J)<=Y2(I) .AND. D1(J)==1 .AND. D2(J)==1 .AND. Y2(J)<TIM(U)) THEN
                  DO K=1,N,1
                     IF (D1(K)==1 .AND. Y1(K)<Y2(J)) THEN
					    VPBRF3(I,U,2)=VPBRF3(I,U,2)*ST(U)*STYG1(K)*RR(J)**2/REAL(N**2,8)
						VPBRF4(I,U,2)=VPBRF4(I,U,2)*ST(U)*STYG2(K)*RR(J)**2/REAL(N**2,8)
					 END IF
			      END DO
			   END IF
			END DO
		 END IF

		 IF (D1(I)==0 .AND. D2(I)==0) THEN
            VPBRF3(I,U,3)=VPBRF3(I,U,3)+ST(U)*SUM(STYG1,D1==1 .AND. Y1(I)<Y1 .AND. Y1<=TIM(U))/REAL(N,8)/W1(I)
		 END IF

		 IF (D2(I)==0) THEN
            VPBRF4(I,U,3)=VPBRF4(I,U,3)+ST(U)*SUM(STYG2,D1==1 .AND. Y2(I)<Y1 .AND. Y1<=TIM(U))/REAL(N,8)/W2(I)
		 END IF


		 DO J=1,N,1
            IF (Y1(I)>=Y1(J) .AND. D1(J)==0 .AND. D2(J)==0) THEN
               DO K=1,N,1
                  IF (D1(K)==1 .AND. Y1(J)<Y1(K) .AND. Y1(K)<=TIM(U)) THEN
                     VPBRF3(I,U,3)=VPBRF3(I,U,3)-ST(U)*STYG1(K)/W1(J)**2/REAL(N**2,8)
				  END IF
			   END DO
			END IF
		 END DO

		 DO J=1,N,1
            IF (Y2(I)>=Y2(J) .AND. D2(J)==0) THEN
               DO K=1,N,1
                  IF (D1(K)==1 .AND. Y2(J)<Y1(K) .AND. Y1(K)<=TIM(U)) THEN
                     VPBRF4(I,U,3)=VPBRF4(I,U,3)-ST(U)*STYG2(K)/W2(J)**2/REAL(N**2,8)
				  END IF
			   END DO
			END IF
		 END DO
	  END DO
      V3(U,1)=SUM((VPBRF3(:,U,1)+VPBRF3(:,U,2)+VPBRF3(:,U,3))**2)/REAL(N**2,8)
      V4(U,1)=SUM((VPBRF4(:,U,1)+VPBRF4(:,U,2)+VPBRF4(:,U,3))**2)/REAL(N**2,8)
      V3(U,2)=SUM(VPBRF3(:,U,1)**2)/REAL(N**2,8)+SUM((VPBRF3(:,U,2)+VPBRF3(:,U,3))**2)/REAL(N**2,8)
	  V4(U,2)=SUM(VPBRF4(:,U,1)**2)/REAL(N**2,8)+SUM((VPBRF4(:,U,2)+VPBRF4(:,U,3))**2)/REAL(N**2,8)
      V3(U,3)=SUM(VPBRF3(:,U,1)**2)/REAL(N**2,8)+SUM(VPBRF3(:,U,2)**2)/REAL(N**2,8)+SUM(VPBRF3(:,U,3)**2)/REAL(N**2,8)
      V4(U,3)=SUM(VPBRF4(:,U,1)**2)/REAL(N**2,8)+SUM(VPBRF4(:,U,2)**2)/REAL(N**2,8)+SUM(VPBRF4(:,U,3)**2)/REAL(N**2,8)
   END DO
   RETURN
END SUBROUTINE YPBF

SUBROUTINE YPBFF(N,Y1,Y2,D1,D2,NT,TIM,PBRF,V1,V2,V3,V4,V5)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'ypbff_' ::YPBFF
   IMPLICIT NONE
   INTEGER,INTENT(IN)::N,NT,D1(N),D2(N)
   REAL(8),INTENT(IN)::Y1(N),Y2(N),TIM(NT)
   REAL(8),INTENT(OUT)::PBRF(NT,5),V1(NT,3),V2(NT,3),V3(NT,3),V4(NT,3),V5(NT,3)

   REAL(8)::ST1(N),ST2(N),G11(N),G12(N),G21(N),G22(N),R1(N),R2(N),STG11(N),STG12(N),STG21(N),STG22(N)
   REAL(8)::ST(NT),G1T(NT),G2T(NT),R(NT),STG1T(NT),STG2T(NT)

   REAL(8),ALLOCATABLE::Y11(:),Y12(:)
   INTEGER,ALLOCATABLE::DD(:)
   INTEGER::I,J,K,U,N1

   REAL(8)::VPBRF1(N,NT,3),VPBRF2(N,NT,3),VPBRF3(N,NT,3),VPBRF4(N,NT,3),VPBRF5(N,NT,3),RR(N)
   REAL(8)::W1(N),W2(N),A1TEMP(NT),A2TEMP(NT),SMIN(NT),S2(NT)
   REAL(8),DIMENSION(N)::A11TEMP,A12TEMP,A21TEMP,A22TEMP


   CALL KM(N,Y1,Y1,1-D1,N,Y1,0,1,G11)
   CALL KM(N,Y1,Y1,1-D1,N,Y2,0,1,G12)
   CALL KM(N,Y2,Y2,1-D2,N,Y1,0,1,G21)
   CALL KM(N,Y2,Y2,1-D2,N,Y2,0,1,G22)
   CALL KM(N,Y1,Y1,1-D1,NT,TIM,0,1,G1T)
   CALL KM(N,Y2,Y2,1-D2,NT,TIM,0,1,G2T)

   CALL KM(N,Y1,Y1,(1-(1-D1)*(1-D2)),NT,TIM,0,1,SMIN)
   CALL KM(N,Y2,Y2,D2,NT,TIM,0,1,S2)

   N1=SUM(D1)
   ALLOCATE(Y11(N1),Y12(N1),DD(N1))
   Y11=PACK(Y1,D1==1)
   Y12=PACK(Y2,D1==1)
   DD=PACK(D2,D1==1)
   CALL KM(N1,Y11,Y12,DD,N,Y1,1,0,ST1)
   CALL KM(N1,Y11,Y12,DD,N,Y2,1,0,ST2)
   CALL KM(N1,Y11,Y12,DD,NT,TIM,1,0,ST)
   A11TEMP=0.0;A12TEMP=0.0;A21TEMP=0.0;A22TEMP=0.0
   DO I=1,N,1
      R2(I)=REAL(COUNT(Y11<=Y2(I) .AND. Y2(I)<=Y12),8)/REAL(N,8)
	  R1(I)=REAL(COUNT(Y11<=Y1(I) .AND. Y1(I)<=Y12),8)/REAL(N,8)
	  W1(I)=REAL(COUNT(Y1>=Y1(I)),8)/REAL(N,8)
	  W2(I)=REAL(COUNT(Y2>=Y2(I)),8)/REAL(N,8)
	  IF (ST1(I)>0.0 .AND. G11(I)>0.0) THEN
	     A11TEMP(I)=R1(I)/(ST1(I)*G11(I))
	  END IF

	  IF (ST2(I)>0.0 .AND. G22(I)>0.0) THEN
	     A22TEMP(I)=R2(I)/(ST2(I)*G22(I))
	  END IF
   END DO
   A1TEMP=0.0;A2TEMP=0.0
   DO I=1,NT,1
      R(I)=REAL(COUNT(Y11<=TIM(I) .AND. TIM(I)<=Y12),8)/REAL(N,8)
      IF (ST(I)>0.0 .AND. G1T(I)>0.0) THEN
	     A1TEMP(I)=R(I)/(ST(I)*G1T(I))
	  END IF

      IF (ST(I)>0.0 .AND. G2T(I)>0.0) THEN
	     A2TEMP(I)=R(I)/(ST(I)*G2T(I))
	  END IF
   END DO

   DEALLOCATE(Y11,Y12,DD)

   PBRF=0.0;V1=1.0;V2=1.0;V3=1.0;V4=1.0;V5=1.0
   VPBRF1=0.0;VPBRF2=0.0;VPBRF3=0.0;VPBRF4=0.0;VPBRF5=0.0
   DO U=1,NT,1
      PBRF(U,5)=S2(U)-SMIN(U)
      IF (G1T(U)>0.0) THEN
         PBRF(U,1)=R(U)/G1T(U)
	  END IF
      IF (G2T(U)>0.0) THEN
         PBRF(U,2)=R(U)/G2T(U)
	  END IF
      PBRF(U,3)=SUM(ST(U)/ST1/G11,ST1>0.0 .AND. G11>0.0 .AND. D1==1 .AND. Y1<=TIM(U))/REAL(N,8)
      PBRF(U,4)=SUM(ST(U)/ST1/G21,ST1>0.0 .AND. G21>0.0 .AND. D1==1 .AND. Y1<=TIM(U))/REAL(N,8)
   END DO
   DO U=1,NT,1
      VPBRF1(:,U,1)=-PBRF(U,1)
      VPBRF2(:,U,1)=-PBRF(U,2)
	  VPBRF3(:,U,1)=-PBRF(U,3)
      VPBRF4(:,U,1)=-PBRF(U,4)

      DO I=1,N,1
         IF (D1(I)==1 .AND. Y1(I)<=TIM(U) .AND. TIM(U)<=Y2(I) .AND. G1T(U)>0.0) THEN
		    VPBRF1(I,U,1)=VPBRF1(I,U,1)+1.0/G1T(U)
		 END IF

		 IF (D1(I)==1 .AND. Y1(I)<=TIM(U) .AND. TIM(U)<=Y2(I) .AND. G2T(U)>0.0) THEN
		    VPBRF2(I,U,1)=VPBRF2(I,U,1)+1.0/G2T(U)
		 END IF

	     IF (D1(I)==1 .AND. Y1(I)<=TIM(U) .AND. ST1(I)>0.0 .AND. G11(I)>0.0) THEN
	        VPBRF3(I,U,1)=VPBRF3(I,U,1)+ST(U)/ST1(I)/G11(I)
         END IF
		 IF (D1(I)==1 .AND. Y1(I)<=TIM(U) .AND. ST1(I)>0.0 .AND. G21(I)>0.0) THEN
            VPBRF4(I,U,1)=VPBRF4(I,U,1)+ST(U)/ST1(I)/G21(I)
		 END IF

		 IF (D1(I)==1 .AND. D2(I)==1 .AND. Y2(I)<=TIM(U) .AND. ST2(I)>0.0 .AND. G12(I)>0.0) THEN
		    VPBRF3(I,U,2)=VPBRF3(I,U,2)-ST(U)/ST2(I)/G12(I)
		 END IF
         IF (D1(I)==1 .AND. D2(I)==1 .AND. Y2(I)<=TIM(U) .AND. ST2(I)>0.0 .AND. G22(I)>0.0) THEN
			VPBRF4(I,U,2)=VPBRF4(I,U,2)-ST(U)/ST2(I)/G22(I)
         END IF

		 IF (D1(I)==1) THEN
            VPBRF3(I,U,2)=VPBRF3(I,U,2)+ST(U)*SUM(1.0/(ST2*G12*R2),ST2>0.0 .AND. G12>0.0 &
            .AND. R2>0.0 .AND. D1==1 .AND. D2==1 .AND. &
            Y1(I)<=Y2 .AND. Y2<=Y2(I) .AND. Y2<=TIM(U))/REAL(N,8)
		 END IF

		 IF (D1(I)==1) THEN
            VPBRF4(I,U,2)=VPBRF4(I,U,2)+ST(U)*SUM(1.0/(ST2*G22*R2),ST2>0.0 .AND. G22>0.0 &
            .AND. R2>0.0 .AND. D1==1 .AND. D2==1 .AND. Y1(I)<=Y2 .AND. Y2<=Y2(I) .AND. Y2<=TIM(U))/REAL(N,8)
		 END IF

		 IF (D1(I)==0 .AND. D2(I)==0 .AND. W1(I)>0.0 .AND. Y1(I)<=TIM(U)) THEN
            VPBRF1(I,U,3)=VPBRF1(I,U,3)+PBRF(U,1)/W1(I)
			VPBRF3(I,U,3)=VPBRF3(I,U,3)+ST(U)*(A1TEMP(U)-A11TEMP(I))/W1(I)
		 END IF

		 IF (D2(I)==0 .AND. W2(I)>0.0 .AND. Y2(I)<=TIM(U)) THEN
		    VPBRF2(I,U,3)=VPBRF2(I,U,3)+PBRF(U,2)/W2(I)
            VPBRF4(I,U,3)=VPBRF4(I,U,3)+ST(U)*(A2TEMP(U)-A22TEMP(I))/W2(I)
		 END IF

         VPBRF1(I,U,3)=VPBRF1(I,U,3)-PBRF(U,1)*SUM(1.0/W1**2,W1>0.0 .AND. D1==0 &
                     .AND. D2==0 .AND. Y1<=Y1(I) .AND. Y1<=TIM(U))/REAL(N,8)
		     VPBRF2(I,U,3)=VPBRF2(I,U,3)-PBRF(U,2)*SUM(1.0/W2**2,W2>0.0 .AND. D2==0 .AND. Y2<=Y2(I) .AND. Y2<=TIM(U))/REAL(N,8)
         VPBRF3(I,U,3)=VPBRF3(I,U,3)-ST(U)*SUM((A1TEMP(U)-A11TEMP)/W1**2,W1>0.0 &
                     .AND. D1==0 .AND. D2==0 .AND. Y1<=Y1(I) .AND. Y1<=TIM(U))/REAL(N,8)
         VPBRF4(I,U,3)=VPBRF4(I,U,3)-ST(U)*SUM((A2TEMP(U)-A22TEMP)/W2**2,W2>0.0 &
                     .AND. D2==0 .AND. Y2<=Y2(I) .AND. Y2<=TIM(U))/REAL(N,8)

         IF (D2(I)==1 .AND. Y2(I)<TIM(U) .AND. W2(I)>0.0) THEN
		    VPBRF5(I,U,1)=VPBRF5(I,U,1)-S2(U)/W2(I)
		 END IF
         VPBRF5(I,U,1)=VPBRF5(I,U,1)+S2(U)*SUM(1.0/W2**2,W2>0.0 .AND. D2==1 .AND. Y2<=Y2(I) .AND. Y2<TIM(U))/REAL(N,8)

         IF ((D1(I)==1 .OR. D2(I)==1) .AND. Y1(I)<TIM(U) .AND. W1(I)>0.0) THEN
		    VPBRF5(I,U,3)=VPBRF5(I,U,3)+SMIN(U)/W1(I)
		 END IF
		 VPBRF5(I,U,3)=VPBRF5(I,U,3)-SMIN(U)*SUM(1.0/W1**2,W1>0.0 .AND. (D1==1 .OR. D2==1) .AND. Y1<=Y1(I) .AND. Y1<TIM(U))/REAL(N,8)

	  END DO
      V3(U,1)=SUM((VPBRF3(:,U,1)+VPBRF3(:,U,2)+VPBRF3(:,U,3))**2)/REAL(N**2,8)
      V4(U,1)=SUM((VPBRF4(:,U,1)+VPBRF4(:,U,2)+VPBRF4(:,U,3))**2)/REAL(N**2,8)
      V3(U,2)=SUM(VPBRF3(:,U,1)**2)/REAL(N**2,8)+SUM((VPBRF3(:,U,2)+VPBRF3(:,U,3))**2)/REAL(N**2,8)
	  V4(U,2)=SUM(VPBRF4(:,U,1)**2)/REAL(N**2,8)+SUM((VPBRF4(:,U,2)+VPBRF4(:,U,3))**2)/REAL(N**2,8)
      V3(U,3)=SUM(VPBRF3(:,U,1)**2)/REAL(N**2,8)+SUM(VPBRF3(:,U,2)**2)/REAL(N**2,8)+SUM(VPBRF3(:,U,3)**2)/REAL(N**2,8)
      V4(U,3)=SUM(VPBRF4(:,U,1)**2)/REAL(N**2,8)+SUM(VPBRF4(:,U,2)**2)/REAL(N**2,8)+SUM(VPBRF4(:,U,3)**2)/REAL(N**2,8)

      V1(U,1)=SUM((VPBRF1(:,U,1)+VPBRF1(:,U,2)+VPBRF1(:,U,3))**2)/REAL(N**2,8)
      V2(U,1)=SUM((VPBRF2(:,U,1)+VPBRF2(:,U,2)+VPBRF2(:,U,3))**2)/REAL(N**2,8)
      V1(U,2)=SUM(VPBRF1(:,U,1)**2)/REAL(N**2,8)+SUM((VPBRF1(:,U,2)+VPBRF1(:,U,3))**2)/REAL(N**2,8)
	  V2(U,2)=SUM(VPBRF2(:,U,1)**2)/REAL(N**2,8)+SUM((VPBRF2(:,U,2)+VPBRF2(:,U,3))**2)/REAL(N**2,8)
      V1(U,3)=SUM(VPBRF1(:,U,1)**2)/REAL(N**2,8)+SUM(VPBRF1(:,U,2)**2)/REAL(N**2,8)+SUM(VPBRF1(:,U,3)**2)/REAL(N**2,8)
      V2(U,3)=SUM(VPBRF2(:,U,1)**2)/REAL(N**2,8)+SUM(VPBRF2(:,U,2)**2)/REAL(N**2,8)+SUM(VPBRF2(:,U,3)**2)/REAL(N**2,8)

      V5(U,1)=SUM((VPBRF5(:,U,1)+VPBRF5(:,U,2)+VPBRF5(:,U,3))**2)/REAL(N**2,8)
      V5(U,2)=SUM(VPBRF5(:,U,1)**2)/REAL(N**2,8)+SUM((VPBRF5(:,U,2)+VPBRF5(:,U,3))**2)/REAL(N**2,8)
      V5(U,3)=SUM(VPBRF5(:,U,1)**2)/REAL(N**2,8)+SUM(VPBRF5(:,U,2)**2)/REAL(N**2,8)+SUM(VPBRF5(:,U,3)**2)/REAL(N**2,8)
   END DO
   RETURN
END SUBROUTINE YPBFF

