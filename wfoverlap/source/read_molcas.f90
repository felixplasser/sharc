!******************************************
!
!    SHARC Program Suite
!
!    Copyright (c) 2019 University of Vienna
!
!    This file is part of SHARC.
!
!    SHARC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHARC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
!
!******************************************

MODULE read_molcas
  USE sysparam
  USE my_alloc
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_molcas_overlap

!--------------------------------------------------------------------------------------------------------------------------------==
!--------------------------------------------------------------------------------------------------------------------------------==
CONTAINS

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_molcas_overlap(file,ovl,nl,nr,same_aos)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: file
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ovl
    INTEGER(KIND=ilong), INTENT(OUT) :: nl
    INTEGER(KIND=ilong), INTENT(OUT) :: nr
    LOGICAL, INTENT(IN) :: same_aos
    INTEGER(KIND=ilong):: nao
    INTEGER(KIND=ilong):: nsym
    INTEGER(KIND=ilong):: irc
    INTEGER(KIND=ilong):: nvalid
    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong):: j
    INTEGER(KIND=ilong):: k
    INTEGER*8:: iCode=6
    INTEGER*8:: ilOne=1
    INTEGER*8:: mc1
    INTEGER(KIND=ilong), DIMENSION(:),ALLOCATABLE :: nbas
    REAL*8, DIMENSION(:),ALLOCATABLE :: temp
    INTEGER:: allocstat

    ! CALL link_it()
    CALL getenvinit()
    CALL fioinit()
    CALL NameRun('RUNFILE')
    CALL Get_iScalar('nSym',nsym)
    IF(nsym.GT.iOne)THEN
      WRITE(0,'("Found ",I2," symmetry groups - this is not a C1-symmetry calculation")'),nsym
      WRITE(6,'("Found ",I2," symmetry groups - this is not a C1-symmetry calculation")'),nsym
      STOP 1
    END IF
    ALLOCATE(nbas(1:nsym))
    CALL Get_iArray('nBas',nbas,nSym)            
    nao = nbas(1)
    
    ! Determine which part of the overlap matrix should be read
    IF(same_aos)THEN
      nl = nao
      nr = nao
    ELSE
      IF ((nl.EQ.-1).AND.(nr.EQ.-1)) THEN
        ! If neither nao_a nor nao_b were specified, assume that they are the same
        nl = nao / 2
        nr = nao / 2
      ELSE
        IF(nl+nr.NE.nao)THEN
          WRITE(6,'("Inconsistent number of AOs given: ", 3I10)') nl, nr, nao
          WRITE(6,*) "file: ",trim(adjustl(file))
          WRITE(0,'("Inconsistent number of AOs given: ", 3I10)') nl, nr, nao
          WRITE(0,*) "file: ",trim(adjustl(file))
          STOP 401
        END IF
      END IF
    END IF

    allocstat=myalloc(ovl,nl,nr,'+ ovl mat')
    IF(allocstat.NE.0) STOP

    allocstat=myalloc(temp,nao*(nao+1)/2)
    IF(allocstat.NE.0) STOP

    mc1=1
    CALL rdonexx(irc,iCode,'Mltpl  0',mc1,temp,ilOne,nvalid,trim(adjustl(file)))

    k=0
    IF (same_aos)THEN
      DO i=1,nbas(1)
        DO j=1,i
          k=k+1
          ovl(i,j)=temp(k)
          ovl(j,i)=temp(k)
        END DO ! j = 1,i
      END DO ! i = 1,nbas
    ELSE
      ! Read only the mixed block
      DO i=1,nbas(1)
        DO j=1,i
          k=k+1
          IF((i > nl).AND.(j <= nl)) ovl(j,i-nl)=temp(k)
        END DO ! j = 1,i
      END DO ! i = 1,nbas
    END IF

    IF(debug)THEN
      WRITE(6,*) "Overlap matrix read from Molcas ",trim(adjustl(file))
      WRITE(6,*)nl, nr
      DO i=1,nl
        DO j=1,nr
          WRITE(6,'(x,ES24.14E3)',ADVANCE='no')ovl(i,j)
        END DO
        WRITE(6,*)
      END DO
    END IF ! debug
    
    allocstat = mydealloc(temp)
    IF(allocstat.NE.0)THEN
      WRITE(6,*) "Problem when deallocating temp matrix in read_mc_overlap - ignoring"
      WRITE(0,*) "Problem when deallocating temp matrix in read_mc_overlap - ignoring"
    END IF

  END SUBROUTINE read_molcas_overlap


! MOLCAS WRAPPER FOR READING AO INTEGRAL FILES TAKEN FROM COLUMBUS

  SUBROUTINE rdonexx(i_rc,i_Option,InLab,i_Comp,Data,i_SymLab,ret_data,fname)
!
!*********************************************************************
!                                                                      *
!   (c) Copyright 1993 by the authors of MOLCAS. All rights reserved   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Purpose: Read data from one-electron integral file               *
!                                                                      *
!     Calling parameters:                                              *
!     rc      : Return code                                            *
!     Option  : Switch to set options                                  *
!     InLab   : Identifier for the data to read                        *
!     Comp    : Composite identifier to select components              *
!     Data    : contains on output the requested data                  *
!     SymLab  : symmetry label of the requested data                   *
!     ret_data: valid data in field Data (DP units)
!                                                                      *
!     Global data declarations (Include file):                         *
!     Parm    : Table of paramaters                                    *
!     rcParm  : Table of return codes                                  *
!     Switch  : Table of options                                       *
!    Common  : Common block containing ToC                            *
!     Data    : Data definitions                                       *
!                                                                      *
!     Local data declarations:                                         *
!     Label   : character*8, used to covert incoming names             *
!     TmpBuf  : I/O buffer                                             *
!     HldBuf  : I/O buffer                                             *
!                                                                     *
!---------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark and M.P. Fuelscher                                  *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!*********************************************************************
!@ifdef molcas_int64
       Implicit Integer*8 (A-Z)
!@ifndef int64
!        integer*4 i_rc,i_option,i_comp,i_listlabels,i_symlab
!@endif
!@else
!      Implicit Integer (A-Z)
!@endif
!----------------------------------------------------------------------*
!                                                                      *
! Return codes:                                                        *
!   rc0000 - No error                                                  *
!   rcOP01 - file is already opened                                    *
!   rcOP02 - file specified as old is not existent                     *
!   rcOP03 - invalid file identifier                                   *
!   rcOP04 - unknown option has been specified                         *
!   rcCL01 - file is not opened                                        *
!   rcRD01 - file is not opened                                        *
!   rcRD02 - illegal options                                           *
!   rcRD03 - information not available                                 *
!   rcRD04 - nSym not defined                                          *
!   rcWR01 - file is not opened                                        *
!   rcWR02 - nSym<1 or nSym>8                                          *
!   rcWR03 - nSym note defined                                         *
!   rcWR04 - Sum(nBas(iSym))>2*mxBas                                   *
!   rcWR05 - Sum(nBas(iSym))<1                                         *
!   rcWR06 - Min(nBas(iSym)<0                                          *
!   rcWR07 - Max(nBas(iSym)>mxBas                                      *
!   rcWR08 - nAtm<1 or nAtm>mxAtm                                      *
!   rcWR09 - nAtm not defined                                          *
!   rcWR10 - nBas not defined                                          *
!   rcWR11 - to many operators                                         *
!                                                                      *
!----------------------------------------------------------------------*
      Parameter (rc0000 = 0)
      Parameter (rcOP01 = rc0000+1)
      Parameter (rcOP02 = rcOP01+1)
      Parameter (rcOP03 = rcOP02+1)
      Parameter (rcOP04 = rcOP03+2)
      Parameter (rcCL01 = rcOP04+1)
      Parameter (rcRD01 = rcCL01+1)
      Parameter (rcRD02 = rcRD01+1)
      Parameter (rcRD03 = rcRD02+1)
      Parameter (rcRD04 = rcRD03+1)
      Parameter (rcWR01 = rcRD04+1)
      Parameter (rcWR02 = rcWR01+1)
      Parameter (rcWR03 = rcWR02+1)
      Parameter (rcWR04 = rcWR03+1)
      Parameter (rcWR05 = rcWR04+1)
      Parameter (rcWR06 = rcWR05+1)
      Parameter (rcWR07 = rcWR06+1)
      Parameter (rcWR08 = rcWR07+1)
      Parameter (rcWR09 = rcWR08+1)
      Parameter (rcWR10 = rcWR09+1)
      Parameter (rcWR11 = rcWR10+1)
!
!ifordeck OneFlags.inc $Revision: 2.26.8.4 $
!----------------------------------------------------------------------*
!                                                                      *
! Switches for one electron integral handlers:                         *
!   sOpSiz - Return only size of array                                 *
!   sNoOri - Do not read origin of operator                            *
!   sNoNuc - Do not read nuclear contribution                          *
!   sRdFst - Read first operator                                       *
!   sRdNxt - Read next operator                                        *
!   sNew   - New file, create.                                         *
!   sXXXXX -                                                           *
!                                                                      *
!----------------------------------------------------------------------*
      Parameter (sOpSiz = 1)
      Parameter (sNoOri = 2*sOpSiz)
      Parameter (sNoNuc = 2*sNoOri)
      Parameter (sRdFst = 2*sNoNuc)
      Parameter (sRdNxt = 2*sRdFst)
      Parameter (sRdCur = 2*sRdNxt)
      Parameter (sXXXXX = 2*sRdCur)
      Parameter (sNew   = 1)
!comdeck OneDat.INC $Revision: 2.26.8.4 $
!----------------------------------------------------------------------*
!                                                                      *
!     Define I/O options                                               *
!                                                                      *
!----------------------------------------------------------------------*
      Parameter ( sWrite = 1          )
      Parameter ( sRead  = 2          )
!----------------------------------------------------------------------*
!                                                                      *
!     Define data conversion factors (machine dependent)               *
!                                                                      *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
       Integer    ItoB,      RtoB,      RtoI
!@ifdef molcas_int64
       Parameter( ItoB = 8 , RtoB = 8 , RtoI = 1  )
!@else
!      Parameter( ItoB = 4 , RtoB = 8 , RtoI = 2  )
!@endif
!
       Parameter (mxAtom = 500)
       Parameter(maxbfn=5000)
       Parameter (mxroot = 100)
       Parameter (mxact  =  50)
       Parameter (mxina  = 100)
       Parameter (mxbas = maxbfn)
       Parameter (mxOrb = maxbfn)
       Parameter (mxSym = 8)
!----------------------------------------------------------------------*
!                                                                      *
!     Define Common /OneDat/                                           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Parameters:                                                          *
!   NaN    - Not a Number = Variable undefined                         *
!   NotNaN - A Number = Variable is defined                            *
!   nAuxDt - extra space for origin and nuclear contribution.          *
!   nTitle - length of title.                                          *
!   MxOp   - Max number of operators                                   *
!   LenOp  - Length of TOC field defining operator                     *
!   MxSym  - Max number of symmetry representations                    *
!   MxAtom - Max number of atoms in system                             *
!   MxBas  - Max number of total basis functions                       *
!   PhyRc  - physical buffer size for DAFILE                           *
!   nBuf   - logical internal buffer size for reading/writing matrices *
!                                                                      *
! Pointers:                                                            *
!   pFID   - File identifier                                           *
!   pVersN - Version number                                            *
!   pTitle - Titleof the problem                                       *
!   pOp    - Operator list                                             *
!   pSym   - Number of irred. representations                          *
!   pSymOp - generator of irred. representation                        *
!   pBas   - Number of basis functions per irred. representation       *
!   pAtom  - Atoms in system                                           *
!   pCoord - Coordinates of atoms in system                            *
!   pPot   - Nuclear-Nuclear repulsion                                 *
!   pCoM   - Coordinates of center of mass                             *
!   pCoC   - Coordinates of center of charge                           *
!   pALbl  - Atom labels                                               *
!   pType  - Basis function symmetry labels                            *
!   pChrge - Charge of atoms in system                                 *
!   pOption- Various options - direct, expert mode, ...                *
!   pNext  - Next free record                                          *
!                                                                      *
! Offsets:                                                             *
!   oLabel - Label of operator                                         *
!   oComp  - Component number                                          *
!   oSymLb - Symmetry label of operator                                *
!   oAddr  - Disk address                                              *
!                                                                      *
!----------------------------------------------------------------------*
      Parameter ( NaN=-1        )
      Parameter ( NotNaN=0      )
      Parameter ( nAuxDt=4      )
!    Bug Fix RL 970429
!     Parameter ( nTitle=RtoI*9 )
      Parameter ( nTitle=(72*2)/ItoB)
!
      Parameter ( PhyRec=1024   )
      Parameter ( nBuf=4*PhyRec )
      Parameter ( MxOp=16384    )
      Parameter ( LenOp=5       )
!
      Parameter ( pFID   = 1                       )
      Parameter ( pVersN = pFID   + 1              )
      Parameter ( pTitle = pVersN + 1              )
      Parameter ( pOp    = pTitle + nTitle + 1     )
      Parameter ( pSym   = pOp    + MxOp*LenOp     )
      Parameter ( pSymOp = pSym   + 1              )
      Parameter ( pBas   = pSymOp + MxSym          )
      Parameter ( pAtom  = pBas   + MxSym          )
      Parameter ( pCoord = pAtom  + 1              )
      Parameter ( pPot   = pCoord + 6*MxAtom + 1   )
      Parameter ( pCoM   = pPot   + 2 + 1          )
      Parameter ( pCoC   = pCoM   + 6 + 1          )
      Parameter ( pALbl  = pCoC   + 6 + 1          )
      Parameter ( pType  = pALbl  + MxAtom + 1     )
      Parameter ( pChrge = pType  + 4*MxBas + 1    )
      Parameter ( pIndex = pChrge + 2*MxAtom + 1   )
      Parameter ( pNext  = pIndex + 2*MxAtom + 1   )
      Parameter ( pOption= pNext  + 1              )
      Parameter ( pEnd   = pOption+ 1              )
      Parameter ( lToc   = 1024*((pEnd+1022)/1024) )
!
      Parameter ( oLabel = 0          )
      Parameter ( oComp  = oLabel + 2 )
      Parameter ( oSymLb = oComp  + 1 )
      Parameter ( oAddr  = oSymLb + 1 )
!
      Parameter ( pLu    = 1          )
      Parameter ( pOpen  = pLu    + 1 )
      Parameter ( lAux   = pOpen  + 1 )
!
      Dimension AuxOne(lAux)
      Dimension TocOne(lToc)
!@ifdef molcas_ext
!      Common /OneDat_/ AuxOne,TocOne
!@else
      Common /OneDat/ AuxOne,TocOne
!@endif
!
      Character*(*) InLab
      !Dimension Data(*)
      REAL*8, DIMENSION(:),ALLOCATABLE :: Data !only required change to original routine AEAF
      Dimension nBas(8)
      CHARACTER*(*) fname
!
      Character*8 TmpLab,Label
      Dimension LabTmp(2)
      Equivalence (TmpLab,LabTmp)
!
      Parameter (lBuf=1024)
      Real*8    TmpBuf(lBuf),AuxBuf(4)
      Logical debug, Close
!
      integer*4 i1,i2,i3,i4,i5,ii
!
      Data CurrOp/1/
      Save CurrOp
!----------------------------------------------------------------------*
!     Start procedure:                                                 *
!     Define inline function (symmetry multiplication)                 *
!----------------------------------------------------------------------*
      MulTab(i,j)=iEor(i-1,j-1)+1
!
      rc=i_rc
      Option=i_option
      comp=i_comp
      symlab=i_symlab
      listlabels=i_listlabels
!
!----------------------------------------------------------------------*
!     Pick up the file definitions                                     *
!----------------------------------------------------------------------*
!@ifdef molcas_ext
!      Call qEnter_('RdOne')
!@else
      Call qEnter('RdOne')
!@endif
      rc    = rc0000
      LuOne = AuxOne(pLu  )
      Open  = AuxOne(pOpen)
!----------------------------------------------------------------------*
!     Check the file status                                            *
!----------------------------------------------------------------------*
      Close=.False.
      If ( Open.ne.1 ) Then
!
!------- Well, I'll open and close it for you under the default name
!
         LuOne=77
!@ifdef molcas_ext
!         LuOne=isFreeUnit_(LuOne)
!@else
         LuOne=isFreeUnit(LuOne)
!@endif
         if (len_trim(fname).gt.8) then
            stop 'RdOneXX: fname > 8'
         elseif (len_trim(fname).lt.8) then
            Label=FNAME(1:len_trim(fname))
            do ii=len_trim(fname)+1,8
            label(ii:ii)=' '
            enddo
         else
         Label=FNAME(1:8)
         endif
         iRC=-1
         iOpt=0
!@ifdef molcas_ext
!         Call OpnOne_(iRC,iOpt,Label,LuOne)
!@else
         Call OpnOne(iRC,iOpt,Label,LuOne)
!@endif
         If (iRC.ne.0) Then
            Write (6,*) 'RdOne: Error opening file'
!@ifdef molcas_ext
!            Call Abend_
!@else
            Call Abend
!@endif
         End If
         Close=.True.
      End If                            
!----------------------------------------------------------------------*
!     Truncate the label to 8 characters and convert it to upper case  *
!----------------------------------------------------------------------*
      Label=InLab
!@ifdef molcas_ext
!      Call UpCase_(Label)
!@else
      Call UpCase(Label)
!@endif
      TmpLab=Label
!----------------------------------------------------------------------*
!     Print debugging information                                      *
!----------------------------------------------------------------------*
      debug=.false.
      If(iAnd(int(option,8),int(1024,8)).ne.0) debug=.true.
      If(debug) Then
         Write(6,*) '<<< Entering RdOne >>>'
         Write(6,'(a,z8)') ' rc on entry:     ',rc
         Write(6,'(a,a)')  ' Label on entry:  ',Label
         Write(6,'(a,z8)') ' Comp on entry:   ',Comp
         Write(6,'(a,z8)') ' SymLab on entry: ',SymLab
         Write(6,'(a,z8)') ' Option on entry: ',Option
      End If
!----------------------------------------------------------------------*
!     Check reading mode                                               *
!----------------------------------------------------------------------*
      If((iAnd(iAnd(option,sRdFst),sRdNxt)).ne.0) then
         Write (6,*) 'RdOne: Invalid option(s)'
         Write (6,*) 'option=',option
!@ifdef molcas_ext
!         Call QTrace_()
!         Call Abend_()
!@else
         Call QTrace()
         Call Abend()
!@endif
      Else If((iAnd(iAnd(option,sRdFst),sRdCur)).ne.0) then
         Write (6,*) 'RdOne: Invalid option(s)'
         Write (6,*) 'option=',option
!@ifdef molcas_ext
!         Call QTrace_()
!         Call Abend_()
!@else
         Call QTrace()
         Call Abend()
!@endif
      Else If((iAnd(iAnd(option,sRdNxt),sRdCur)).ne.0) then
         Write (6,*) 'RdOne: Invalid option(s)'
         Write (6,*) 'option=',option
!@ifdef molcas_ext
!         Call QTrace_()
!         Call Abend_()
!@else
         Call QTrace()
         Call Abend()
!@endif
      End If
!----------------------------------------------------------------------*
!     Load back TocOne                                                 *
!----------------------------------------------------------------------*
      iDisk=0
      ia1=2
!@ifdef molcas_ext
!      Call iDaFile_(LuOne,ia1,TocOne,lToc,iDisk)
!@else
      Call iDaFile(LuOne,ia1,TocOne,lToc,iDisk)
!@endif
!----------------------------------------------------------------------*
!     Read data from ToC                                               *
!----------------------------------------------------------------------*
      NoGo=sRdFst+sRdNxt+sRdCur
!----------------------------------------------------------------------*
!     Read operators from integral records                             *
!----------------------------------------------------------------------*
      If (iAnd(option,sRdNxt).ne.0) Then
         CurrOp=CurrOp+1
         If (CurrOp.gt.MxOp) Then
            CurrOp=0
         Else If(TocOne(pOp+LenOp*(CurrOp-1)+oLabel).eq.Nan) Then
            CurrOp=0
         Else
            i=CurrOp
            LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
!@ifndef molcas_int64
!c           due to equivalence statement between character*8 and integer*4
!            LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1 )
!@endif
            Label=TmpLab
            InLab=Label
            SymLab=TocOne(pOp+LenOp*(i-1)+oSymLb)
            Comp=TocOne(pOp+LenOp*(i-1)+oComp   )
         End If
      Else If(iAnd(option,sRdFst).ne.0) Then
         CurrOp=1
         If(TocOne(pOp+LenOp*(CurrOp-1)+oLabel).eq.Nan) Then
            CurrOp=0
         Else
            i=CurrOp
            LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
!@ifndef molcas_int64
!            LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1 )
!@endif
            Label=TmpLab
            InLab=Label
            SymLab=TocOne(pOp+LenOp*(i-1)+oSymLb)
            Comp=TocOne(pOp+LenOp*(i-1)+oComp   )
         End If
      Else If(iAnd(option,sRdCur).ne.0) Then
         If(CurrOp.lt.1 .or. CurrOp.gt.MxOp) Then
            CurrOp=0
         Else If(TocOne(pOp+LenOp*(CurrOp-1)+oLabel).eq.Nan) Then
            CurrOp=0
         Else
            i=CurrOp
            LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
            Label=TmpLab
            InLab=Label
            SymLab=TocOne(pOp+LenOp*(i-1)+oSymLb)
            Comp=TocOne(pOp+LenOp*(i-1)+oComp   )
         End If
      Else
         CurrOp=0
         Do 500 i=MxOp,1,-1
            LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
!@ifndef molcas_int64
!            LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1 )
!@endif
            CmpTmp=TocOne(pOp+LenOp*(i-1)+oComp   )
            TmpCmp=Comp
            If(TmpLab.eq.Label .and. CmpTmp.eq.TmpCmp) CurrOp=i
            if (debug) then 
               if (CurrOp.ne.0) then
               if (debug) then 
                write(6,*) 'picked Label',Tmplab,' Component',cmptmp,'currop=',currop
               endif
                goto 501
               endif     
            endif
  500      Continue
  501      continue
      End If
      If(CurrOp.eq.0) Then
         rc=rcRD03
         if (debug) then
         Write (*,*) 'RdOne: Information not available'
         Write (*,*) 'Option=',Option
         Write (*,*) 'Comp=',Comp
         Write (*,*) 'SymLab=',SymLab
         Write (*,*) 'Label=',Label
         endif 
         Go To 999
      End If
      SymLab=TocOne(pOp+LenOp*(CurrOp-1)+oSymLb)
      If(debug) Then
        write(6,*) 'Computed Symlab=',Symlab
      endif
      Len=0
!@ifdef molcas_ext
!      Call Get_iScalar_('nSym',in)
!      Call Get_iArray_('nBas',nBas,in)
!@else
      Call Get_iScalar('nSym',in)
      Call Get_iArray('nBas',nBas,in)
!@endif
      Do 510 i=1,in
      Do 510 j=1,i
         ij=MulTab(i,j)-1
         If(iAnd(2**ij,SymLab).ne.0) Then
            If(i.eq.j) Then
               Len=Len+nBas(i)*(nBas(i)+1)/2
            Else
               Len=Len+nBas(i)*nBas(j)
            End If
         End If
  510   Continue
      Data(1)=Len
      If ( IAND(option,sOpSiz).eq.0 ) Then
         IndAux = 0
         IndDta = 0
         iDisk=TocOne(pOp+LenOp*(CurrOp-1)+oAddr)
         Do i = 0,Len+3,lBuf
           nCopy  = MAX(int(0,8),MIN(lBuf,Len+int(4,8)-i))
           nSave  = MAX(int(0,8),MIN(lBuf,Len-i))
           ia1=2
!@ifdef molcas_ext
!           Call dDaFile_(LuOne,ia1,TmpBuf,nCopy,iDisk)
!@else
           Call dDaFile(LuOne,ia1,TmpBuf,nCopy,iDisk)
!@endif
! Problem: molcas is compiled with fortran (integer*8) Blas
!          or blas wrapper,
!          Columbus is here using essl with integer*4 args
!          and links in the "wrong library"
!@if defined ( molcas_int64 ) && ( ! defined ( int64 ))
!           i1=nSave
!           i2=1
!           call dcopy_wr(i1,tmpBuf,i2,Data(IndDta+1),i2)
!@else
           Call dCopy_wr(nSave,TmpBuf,1,Data(IndDta+1),1)
!@endif
           IndDta = IndDta+RtoI*nSave
           Do j = nSave+1,nCopy
             IndAux = IndAux+1
            AuxBuf(IndAux) = TmpBuf(nSave+IndAux)
             AuxBuf(IndAux) = TmpBuf(j)
           End Do
         End Do
         ret_data=inddta 
         If(iAnd(sNoOri,option).eq.0) Then
!@if defined (molcas_int64) && (! defined (int64))
!           i1=3
!           i2=1
!           call dcopy_wr(i1,Auxbuf,i2,Data(IndDta+1),i2)
!@else
           Call dCopy_wr(3,AuxBuf,1,Data(IndDta+1),1)
           if (debug) then
             write(6,'(a,3f12.6)') 'Orign of Operator:',(Auxbuf(i1),i1=1,3)
           endif
!@endif
         ret_data=inddta+RtoI*3 
         End If
         If(iAnd(sNoNuc,option).eq.0) Then
!@if defined (molcas_int64) && (! defined (int64))
!            i1=1 
!            i2=1
!            call dcopy_wr(i1,Auxbuf(4),i2,Data(IndDta+RtoI*3+1),i2)
!@else
           Call dCopy_wr(1,AuxBuf(4),1,Data(IndDta+RtoI*3+1),1)
           if (debug) then
             write(6,'(a,f12.6)') 'Nuclear Contribution:',Auxbuf(4)
           endif
         ret_data=inddta+RtoI*3+RtoI*1 
!@endif
         End If
      End If
!
 999  Continue
        ret_data=ret_data/RtoI  
!      in DP units
      If (Close) Then
!       Write (*,*) ' I will close the file for you!'
         iRC=-1
         iOpt=0
!@ifdef molcas_ext
!         Call ClsOne_(iRC,iOpt)
!@else
         Call ClsOne(iRC,iOpt)
!@endif
         If (iRC.ne.0) Then
            Write (6,*) 'RdOne: Error closing file'
!@ifdef molcas_ext
!            Call Abend_
!@else
            Call Abend
!@endif
         End If
      End If              
!----------------------------------------------------------------------*
!     Terminate procedure                                              *
!----------------------------------------------------------------------*
      i_rc=int(rc,4)
      i_symlab=int(symlab,4)
!@ifdef molcas_ext
!      Call qExit_('RdOne')
!@else
      Call qExit('RdOne')
!@endif
      Return
      END SUBROUTINE rdonexx


! DALTON/COLUMBUS ROUTINE REQUIRED FOR rdonexx 

      SUBROUTINE dcopy_wr(N,DX,INCX,DY,INCY)
!
!     COPY DOUBLE PRECISION DX TO DOUBLE PRECISION DY.
!     FOR I = 0 TO N-1, COPY DX(LX+I*INCX) TO DY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.
!
      DOUBLE PRECISION DX(1),DY(1)
      !INTEGER DY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N,NS 
#if defined (VAR_SBLAS)
       call scopy(n,dx,incx,dy,incy)
       return
#endif 
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
!
!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
!
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
!
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
!
   60 CONTINUE
      NS=N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DX(I)
   70     CONTINUE
      RETURN
      END SUBROUTINE dcopy_wr

  !--------------------------------------------------------------------------------------------------------------------------------

END MODULE read_molcas
