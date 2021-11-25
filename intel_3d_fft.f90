!  This file contains:
!  1.  SUBROUTINE fftc2c
!  2.  SUBROUTINE ifftc2c
!  3.  SUBROUTINE fftc2c_inplace
!  4.  SUBROUTINE ifftc2c_inplace
!  5.  SUBROUTINE fftr2c_inplace
!  6.  SUBROUTINE ifftc2r_inplace
!  7.  SUBROUTINE fftr2c
!  8.  SUBROUTINE ifftc2r
!  
!  CREATOR:  Yu Yu @ SHAO 2011.10.29
!  LAST MODIFIER:  Yu Yu @ SJTU 2014.08.31

!======================================
SUBROUTINE fftc2c(indata1d,outdata1d,L)
!  This SUBROUTINE compute forward FFT from COMPLEX to COMPLEX
!  with forward scale factor 1.0
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
COMPLEX(4)::indata1d(L(1)*L(2)*L(3))
COMPLEX(4)::outdata1d(L(1)*L(2)*L(3))
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::status

print*,'begin fftc2c'

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, L )
status = DftiSetValue( my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputeForward( my_desc_handle, indata1d, outdata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'fftc2c success'

END SUBROUTINE fftc2c

!=======================================
SUBROUTINE ifftc2c(indata1d,outdata1d,L)
!  This SUBROUTINE compute backward FFT from COMPLEX to COMPLEX
!  with forward scale factor 1.0/N_grid
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
COMPLEX(4)::indata1d(L(1)*L(2)*L(3))
COMPLEX(4)::outdata1d(L(1)*L(2)*L(3))
REAL(4)::factor
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::status

print*,'begin ifftc2c'

factor=1.0/L(1)/L(2)/L(3)

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, L )
status = DftiSetValue( my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
status = DftiSetValue( my_desc_handle, DFTI_BACKWARD_SCALE, factor )
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputeBackward( my_desc_handle, indata1d, outdata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'ifftc2c success'

END SUBROUTINE ifftc2c

!======================================
SUBROUTINE fftc2c_inplace(indata1d,L)
!  This SUBROUTINE compute forward FFT from COMPLEX to COMPLEX
!  with forward scale factor 1.0
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
COMPLEX(4)::indata1d(L(1)*L(2)*L(3))
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::status

print*,'begin fftc2c'

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, L )
status = DftiSetValue( my_desc_handle, DFTI_PLACEMENT, DFTI_INPLACE )
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputeForward( my_desc_handle, indata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'fftc2c success'

END SUBROUTINE fftc2c_inplace

!=======================================
SUBROUTINE ifftc2c_inplace(indata1d,L)
!  This SUBROUTINE compute backward FFT from COMPLEX to COMPLEX
!  with forward scale factor 1.0/N_grid
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
COMPLEX(4)::indata1d(L(1)*L(2)*L(3))
REAL(4)::factor
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::status

print*,'begin ifftc2c'

factor=1.0/L(1)/L(2)/L(3)

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, L )
status = DftiSetValue( my_desc_handle, DFTI_PLACEMENT, DFTI_INPLACE )
status = DftiSetValue( my_desc_handle, DFTI_BACKWARD_SCALE, factor )
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputeBackward( my_desc_handle, indata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'ifftc2c success'

END SUBROUTINE ifftc2c_inplace 

!=======================================
SUBROUTINE fftr2c_inplace(indata1d,L)
!  This SUBROUTINE compute forward FFT from REAL to COMPLEX
!  with forward scale factor 1.0
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
REAL(4)::indata1d((L(1)+2)*L(2)*L(3))
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::strides_in(4)
INTEGER::strides_out(4)
INTEGER::status

strides_in(1)=0
strides_in(2)=1
strides_in(3)=L(1)+2
strides_in(4)=(L(1)+2)*L(2)
strides_out(1)=0
strides_out(2)=1
strides_out(3)=L(1)/2+1
strides_out(4)=(L(1)/2+1)*L(2)

print*,'begin fftr2c'

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_REAL, 3, L )
Status = DftiSetValue( my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
Status = DftiSetValue( my_desc_handle, DFTI_INPUT_STRIDES, strides_in)
Status = DftiSetValue( my_desc_handle, DFTI_OUTPUT_STRIDES, strides_out)
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputeForward( my_desc_handle, indata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'fftr2c success'

END SUBROUTINE fftr2c_inplace

!=======================================
SUBROUTINE ifftc2r_inplace(indata1d,L)
!  This SUBROUTINE compute forward FFT from REAL to COMPLEX
!  with forward scale factor 1.0
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
REAL(4)::indata1d((L(1)+2)*L(2)*L(3))
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::strides_in(4)
INTEGER::strides_out(4)
INTEGER::status

strides_in(1)=0
strides_in(2)=1
strides_in(3)=L(1)/2+1
strides_in(4)=(L(1)/2+1)*L(2)
strides_out(1)=0
strides_out(2)=1
strides_out(3)=L(1)+2
strides_out(4)=(L(1)+2)*L(2)

print*,'begin ifftc2r'

factor=1.0/L(1)/L(2)/L(3)

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_REAL, 3, L )
Status = DftiSetValue( my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
Status = DftiSetValue( my_desc_handle, DFTI_INPUT_STRIDES, strides_in)
Status = DftiSetValue( my_desc_handle, DFTI_OUTPUT_STRIDES, strides_out)
status = DftiSetValue( my_desc_handle, DFTI_BACKWARD_SCALE, factor )
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputebackward( my_desc_handle, indata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'ifftc2r success'

END SUBROUTINE ifftc2r_inplace

!=======================================
SUBROUTINE fftr2c(indata1d,outdata1d,L)
!  This SUBROUTINE compute forward FFT from REAL to COMPLEX
!  with forward scale factor 1.0
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
REAL(4)::indata1d((1)*L(2)*L(3))
REAL(4)::outdata1d((L(1)+2)*L(2)*L(3))
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::strides_in(4)
INTEGER::strides_out(4)
INTEGER::status

strides_in(1)=0
strides_in(2)=1
strides_in(3)=L(1)
strides_in(4)=L(1)*L(2)
strides_out(1)=0
strides_out(2)=1
strides_out(3)=L(1)/2+1
strides_out(4)=(L(1)/2+1)*L(2)

print*,'begin fftr2c'

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_REAL, 3, L )
status = DftiSetValue( my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
Status = DftiSetValue( my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
Status = DftiSetValue( my_desc_handle, DFTI_INPUT_STRIDES, strides_in)
Status = DftiSetValue( my_desc_handle, DFTI_OUTPUT_STRIDES, strides_out)
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputeForward( my_desc_handle, indata1d, outdata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'fftr2c success'

END SUBROUTINE fftr2c

!=======================================
SUBROUTINE ifftc2r(indata1d,outdata1d,L)
!  This SUBROUTINE compute forward FFT from REAL to COMPLEX
!  with forward scale factor 1.0
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
REAL(4)::indata1d((L(1)+2)*L(2)*L(3))
REAL(4)::outdata1d(L(1)*L(2)*L(3))
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::strides_in(4)
INTEGER::strides_out(4)
INTEGER::status

strides_in(1)=0
strides_in(2)=1
strides_in(3)=L(1)/2+1
strides_in(4)=(L(1)/2+1)*L(2)
strides_out(1)=0
strides_out(2)=1
strides_out(3)=L(1)
strides_out(4)=L(1)*L(2)

print*,'begin ifftc2r'

factor=1.0/L(1)/L(2)/L(3)

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_REAL, 3, L )
status = DftiSetValue( my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
Status = DftiSetValue( my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
Status = DftiSetValue( my_desc_handle, DFTI_INPUT_STRIDES, strides_in)
Status = DftiSetValue( my_desc_handle, DFTI_OUTPUT_STRIDES, strides_out)
status = DftiSetValue( my_desc_handle, DFTI_BACKWARD_SCALE, factor )
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputebackward( my_desc_handle, indata1d, outdata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'ifftc2r success'

END SUBROUTINE ifftc2r






!EOF
