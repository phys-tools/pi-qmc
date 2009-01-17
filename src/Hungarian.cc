// $Id: Hungarian.cc,v 1.3 2006/10/18 17:08:18 jshumwa Exp $
/*  Copyright (C) 2004-2006 John B. Shumway, Jr.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/*      SUBROUTINE ASSCT ( N, A, C, T )                            
      INTEGER A(130,131), C(130), CH(130), LC(130), LR(130),
     *        LZ(130), NZ(130), RH(131), SLC(130), SLR(130),
     *        U(131)
      INTEGER H, Q, R, S, T
      EQUIVALENCE (LZ,RH), (NZ,CH)
*/
// THIS SUBROUTINE SOLVES THE SQUARE ASSIGNMENT PROBLEM
// THE MEANING OF THE INPUT PARAMETERS IS
// N = NUMBER OF ROWS AND COLUMNS OF THE COST MATRIX, WITH
//     THE CURRENT DIMENSIONS THE MAXIMUM VALUE OF N IS 130
// A(I,J) = ELEMENT IN ROW I AND COLUMN J OF THE COST MATRIX
// ( AT THE END OF COMPUTATION THE ELEMENTS OF A ARE CHANGED)
// THE MEANING OF THE OUTPUT PARAMETERS IS
// C(J) = ROW ASSIGNED TO COLUMN J (J=1,N)
// T = COST OF THE OPTIMAL ASSIGNMENT
// ALL PARAMETERS ARE INTEGER
// THE MEANING OF THE LOCAL VARIABLES IS
// A(I,J) = ELEMENT OF THE COST MATRIX IF A(I,J) IS POSITIVE,
//          COLUMN OF THE UNASSIGNED ZERO FOLLOWING IN ROW I
//          (I=1,N) THE UNASSIGNED ZERO OF COLUMN J (J=1,N)
//          IF A(I,J) IS NOT POSITIVE
// A(I,N+1) = COLUMN OF THE FIRST UNASSIGNED ZERO OF ROW I
//            (I=1,N)
// CH(I) = COLUMN OF THE NEXT UNEXPLORED AND UNASSIGNED ZERO
/*C         OF ROW I (I=1,N)
C LC(J) = LABEL OF COLUMN J (J=1,N)
C LR(I) = LABEL OF ROW I (I=1,N)
C LZ(I) = COLUMN OF THE LAST UNASSIGNED ZERO OF ROW I(I=1,N)
C NZ(I) = COLUMN OF THE NEXT UNASSIGNED ZERO OF ROW I(I=1,N)
C RH(I) = UNEXPLORED ROW FOLLOWING THE UNEXPLORED ROW I
C         (I=1,N)
C RH(N+1) = FIRST UNEXPLORED ROW
C SLC(K) = K-TH ELEMENT CONTAINED IN THE SET OF THE LABELLED
C          COLUMNS
C SLR(K) = K-TH ELEMENT CONTAINED IN THE SET OF THE LABELLED
C          ROWS
C U(I) = UNASSIGNED ROW FOLLOWING THE UNASSIGNED ROW I
C        (I=1,N)
C U(N+1) = FIRST UNASSIGNED ROW
C
C THE VECTORS C,CH,LC,LR,LZ,NZ,SLC,SLR MUST BE DIMENSIONED
C AT LEAST AT (N), THE VECTORS RH,U AT  LEAST AT (N+1),
C THE MATRIX A AT LEAST AT (N,N+1)
C
C INITIALIZATION
      MAXNUM = 10**14
      NP1 = N+1
      DO 10 J=1,N
        C(J) = 0
        LZ(J) = 0
        NZ(J) = 0
        U(J) = 0
   10 CONTINUE
      U(NP1) = 0
      T = 0
C REDUCTION OF THE INITIAL COST MATRIX
      DO 40 J=1,N
        S = A(1,J)
        DO 20 L=2,N
          IF ( A(L,J) .LT. S ) S = A(L,J)
   20   CONTINUE
        T = T+S
        DO 30 I=1,N
          A(I,J) = A(I,J)-S
   30   CONTINUE
   40 CONTINUE
      DO 70 I=1,N
        Q = A(I,1)
        DO 50 L=2,N
          IF ( A(I,L) .LT. Q ) Q = A(I,L)
   50   CONTINUE
        T = T+Q
        L = NP1
        DO 60 J=1,N
          A(I,J) = A(I,J)-Q
          IF ( A(I,J) .NE. 0 ) GO TO 60
          A(I,L) = -J
          L = J
   60   CONTINUE
   70 CONTINUE
C CHOICE OF THE INITIAL SOLUTION
      K = NP1
      DO 140 I=1,N
        LJ = NP1
        J = -A(I,NP1)
   80   IF ( C(J) .EQ. 0 ) GO TO 130
        LJ = J
        J = -A(I,J)
        IF ( J .NE. 0 ) GO TO 80
        LJ = NP1
        J = -A(I,NP1)
   90   R = C(J)
        LM = LZ(R)
        M = NZ(R)
  100   IF ( M .EQ. 0 ) GO TO 110
        IF ( C(M) .EQ. 0 ) GO TO 120
        LM = M
        M = -A(R,M)
      GO TO 100
  110   LJ = J
        J = -A(I,J)
        IF ( J .NE. 0 ) GO TO 90
        U(K) = I
        K = I
        GO TO 140
  120   NZ(R) = -A(R,M)
        LZ(R) = J
        A(R,LM) = -J
        A(R,J) = A(R,M)
        A(R,M) = 0
        C(M) = R
  130   C(J) = I
        A(I,LJ) = A(I,J)
        NZ(I) = -A(I,J)
        LZ(I) = LJ
        A(I,J) = 0
  140 CONTINUE
C RESEARCH OF A NEW ASSIGNMENT
  150 IF ( U(NP1) .EQ. 0 ) RETURN
      DO 160 I=1,N
        CH(I) = 0
        LC(I) = 0
        LR(I) = 0
        RH(I) = 0
  160 CONTINUE
      RH(NP1) = -1
      KSLC = 0
      KSLR = 1
      R = U(NP1)
      LR(R) = -1
      SLR(1) = R
      IF ( A(R,NP1) .EQ. 0 ) GO TO 220
  170 L = -A(R,NP1)
      IF ( A(R,L) .EQ. 0 ) GO TO 180
      IF ( RH(R) .NE. 0 ) GO TO 180
      RH(R) = RH(NP1)
      CH(R) = -A(R,L)
      RH(NP1) = R
  180 IF ( LC(L) .EQ. 0 ) GO TO 200
      IF ( RH(R) .EQ. 0 ) GO TO 210
  190 L = CH(R)
      CH(R) = -A(R,L)
      IF ( A(R,L) .NE. 0 ) GO TO 180
      RH(NP1) = RH(R)
      RH(R) = 0
      GO TO 180
  200 LC(L) = R
      IF ( C(L) .EQ. 0 ) GO TO 360
      KSLC = KSLC+1
      SLC(KSLC) = L
      R = C(L)
      LR(R) = L
      KSLR = KSLR+1
      SLR(KSLR) = R
      IF ( A(R,NP1) .NE. 0 ) GO TO 170
  210 CONTINUE
      IF ( RH(NP1) .GT. 0 ) GO TO 350
C REDUCTION OF THE CURRENT COST MATRIX
  220 H = MAXNUM
      DO 240 J=1,N
        IF ( LC(J) .NE. 0 ) GO TO 240
        DO 230 K=1,KSLR
          I = SLR(K)
          IF ( A(I,J) .LT. H ) H = A(I,J)
  230   CONTINUE
  240 CONTINUE
      T = T+H
      DO 290 J=1,N
        IF ( LC(J) .NE. 0 ) GO TO 290
        DO 280 K=1,KSLR
          I = SLR(K)
          A(I,J) = A(I,J)-H
          IF ( A(I,J) .NE. 0 ) GO TO 280
          IF ( RH(I) .NE. 0 ) GO TO 250
          RH(I) = RH(NP1)
          CH(I) = J
          RH(NP1) = I
  250     L = NP1
  260     NL = -A(I,L)
          IF ( NL .EQ. 0 ) GO TO 270
          L = NL
          GO TO 260
  270     A(I,L) = -J
  280   CONTINUE
  290 CONTINUE
      IF ( KSLC .EQ. 0 ) GO TO 350
      DO 340 I=1,N
        IF ( LR(I) .NE. 0 ) GO TO 340
        DO 330 K=1,KSLC
          J = SLC(K)
          IF ( A(I,J) .GT. 0 ) GO TO 320
          L = NP1
  300     NL = - A(I,L)
          IF ( NL .EQ. J ) GO TO 310
          L = NL
          GO TO 300
  310     A(I,L) = A(I,J)
          A(I,J) = H
          GO TO 330
  320     A(I,J) = A(I,J)+H
  330   CONTINUE
  340 CONTINUE
  350 R = RH(NP1)
      GO TO 190
C ASSIGNMENT OF A NEW ROW
  360 C(L) = R
      M = NP1
  370 NM = -A(R,M)
      IF ( NM .EQ. L ) GO TO 380
      M = NM
      GO TO 370
  380 A(R,M) = A(R,L)
      A(R,L) = 0
      IF ( LR(R) .LT. 0 ) GO TO 390
      L = LR(R)
      A(R,L) = A(R,NP1)
      A(R,NP1) = -L
      R = LC(L)
      GO TO 360
  390 U(NP1) = U(R)
      U(R) = 0
      GO TO 150
      END
*/
