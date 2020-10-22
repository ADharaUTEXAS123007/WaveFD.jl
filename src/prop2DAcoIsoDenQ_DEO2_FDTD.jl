mutable struct Prop2DAcoIsoDenQ_DEO2_FDTD
    p::Ptr{Cvoid}
end

function Prop2DAcoIsoDenQ_DEO2_FDTD(;
        nz=0,
        nx=0,
        nsponge=60,
        nbz=24,
        nbx=12,
        dz=5.0,
        dx=5.0,
        dt=0.001,
        nthreads=Sys.CPU_THREADS,
        freesurface=true,
        freqQ=5.0,
        qMin=0.1,
        qInterior=100.0)
    nz,nx,nsponge,nbz,nbx,nthreads = map(x->round(Int,x), (nz,nx,nsponge,nbz,nbx,nthreads))
    dz,dx,dt = map(x->Float32(x), (dz,dx,dt))

    @assert nx > 0
    @assert nz > 0
    @assert nsponge > 0
    @assert nthreads > 0
    @assert nbx > 0
    @assert nbz > 0

    fs = freesurface ? 1 : 0

    p = ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_alloc, libprop2DAcoIsoDenQ_DEO2_FDTD),
        Ptr{Cvoid},
        (Clong, Clong,    Clong, Clong, Clong,   Cfloat, Cfloat, Cfloat, Clong, Clong),
         fs,    nthreads, nx,    nz,    nsponge, dx,     dz,     dt,     nbx,   nbz)

    ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_SetupDtOmegaInvQ, libprop2DAcoIsoDenQ_DEO2_FDTD),
        Cvoid,
        (Ptr{Cvoid}, Cfloat, Cfloat, Cfloat),
         p,          freqQ,  qMin,   qInterior)

    Prop2DAcoIsoDenQ_DEO2_FDTD(p)
end

free(prop::Prop2DAcoIsoDenQ_DEO2_FDTD) = ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_free, libprop2DAcoIsoDenQ_DEO2_FDTD), Cvoid, (Ptr{Cvoid},), prop.p)

function size(prop::Prop2DAcoIsoDenQ_DEO2_FDTD)
    nz = ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_getNz, libprop2DAcoIsoDenQ_DEO2_FDTD), (Clong), (Ptr{Cvoid},), prop.p)
    nx = ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_getNx, libprop2DAcoIsoDenQ_DEO2_FDTD), (Clong), (Ptr{Cvoid},), prop.p)
    (nz,nx)
end

for _f in (:V, :B, :PSpace, :MSpace, :PCur, :POld, :TmpPx, :TmpPz, :DtOmegaInvQ)
    symf = "Prop2DAcoIsoDenQ_DEO2_FDTD_get" * string(_f)
    @eval $(_f)(prop::Prop2DAcoIsoDenQ_DEO2_FDTD) = unsafe_wrap(Array, ccall(($symf, libprop2DAcoIsoDenQ_DEO2_FDTD), Ptr{Float32}, (Ptr{Cvoid},), prop.p), size(prop), own=false)
end

propagateforward!(prop::Prop2DAcoIsoDenQ_DEO2_FDTD) = ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_TimeStep, libprop2DAcoIsoDenQ_DEO2_FDTD), Cvoid, (Ptr{Cvoid},), prop.p)
propagateadjoint!(prop::Prop2DAcoIsoDenQ_DEO2_FDTD) = propagateforward!(prop) # self-adjoint

scale_spatial_derivatives!(prop::Prop2DAcoIsoDenQ_DEO2_FDTD) =
    ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_ScaleSpatialDerivatives, libprop2DAcoIsoDenQ_DEO2_FDTD), Cvoid, (Ptr{Cvoid},), prop.p)

function forwardBornInjection!(prop::Prop2DAcoIsoDenQ_DEO2_FDTD,dmodelv,wavefieldp)
    ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_ForwardBornInjection, libprop2DAcoIsoDenQ_DEO2_FDTD), Cvoid,
        (Ptr{Cvoid},Ptr{Cfloat},Ptr{Cfloat}),
         prop.p,    dmodelv,    wavefieldp)
end

abstract type ImagingCondition end
struct ImagingConditionStandard <: ImagingCondition end
function adjointBornAccumulation!(prop::Prop2DAcoIsoDenQ_DEO2_FDTD,imagingcondition::ImagingConditionStandard,dmodelv,wavefieldp)
    ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_AdjointBornAccumulation, libprop2DAcoIsoDenQ_DEO2_FDTD), Cvoid,
        (Ptr{Cvoid},Ptr{Cfloat},Ptr{Cfloat}),
         prop.p,    dmodelv,    wavefieldp)
end

# const wavefieldseparationlib=normpath(joinpath(Base.source_path(),"../libprop2DAcoIsoDenQ_DEO2_FDTD"))
# struct ImagingConditionWaveFieldSeparation <: ImagingCondition end
# function adjointBornAccumulation!(prop::Prop2DAcoIsoDenQ_DEO2_FDTD,imagingcondition::ImagingConditionWaveFieldSeparation,dmodelv,wavefieldp)
#     ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_AdjointBornAccumulation_wavefieldsep, wavefieldseparationlib), Cvoid,
#         (Ptr{Cvoid},Ptr{Cfloat},Ptr{Cfloat}),
#          prop.p,    dmodelv,    wavefieldp)
# end

# abstract type ImagingCondition end
# struct ImagingConditionStandard <: ImagingCondition end
# function adjointBornAccumulation!(prop::Prop2DAcoIsoDenQ_DEO2_FDTD,
#         imagingcondition::ImagingConditionStandard,dmodelv,wavefieldp)
#     nz,nx = size(prop)
#     _B = B(prop)
#     _V = V(prop)
#     _POld = POld(prop)

#     @views @. dmodelv += 2 * _B * wavefieldp * _POld /_V^3
#     nothing
# end

struct ImagingConditionWaveFieldSeparation <: ImagingCondition end
function adjointBornAccumulation!(prop::Prop2DAcoIsoDenQ_DEO2_FDTD,
    imagingcondition::ImagingConditionWaveFieldSeparation,dmodelv,wavefieldp)
    nz,nx = size(prop)
    nfft = 2 * nz
    nfft2 = div(nfft,2)
    _B = B(prop)
    _V = V(prop)
    _POld = POld(prop)

    tmp_nlfup = zeros(ComplexF32,nfft,nx)
    tmp_adjdn = zeros(ComplexF32,nfft,nx)

    @views tmp_nlfup[1:nz,:] .= wavefieldp
    @views tmp_adjdn[1:nz,:] .= _POld

    FFTW.set_num_threads(Sys.CPU_THREADS)
    plan_fft_forward = plan_fft!(tmp_nlfup, 1;flags=FFTW.MEASURE)
    
    mul!(tmp_nlfup, plan_fft_forward, tmp_nlfup)
    mul!(tmp_adjdn, plan_fft_forward, tmp_adjdn)

    @views tmp_nlfup[1:nfft2+1,:] .= 0
    @views tmp_adjdn[nfft2:end,:] .= 0

    ldiv!(tmp_nlfup, plan_fft_forward, tmp_nlfup)
    ldiv!(tmp_adjdn, plan_fft_forward, tmp_adjdn)
    
    @strided dmodelv .+= 2 .* _B .* real(tmp_nlfup[1:nz,:] .* tmp_adjdn[1:nz,:]) ./ _V.^3
    nothing
end

function show(io::IO, prop::Prop2DAcoIsoDenQ_DEO2_FDTD)
    nx = ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_getNx, libprop2DAcoIsoDenQ_DEO2_FDTD), (Clong), (Ptr{Cvoid},), prop.p)
    nz = ccall((:Prop2DAcoIsoDenQ_DEO2_FDTD_getNz, libprop2DAcoIsoDenQ_DEO2_FDTD), (Clong), (Ptr{Cvoid},), prop.p)
    write(io, "Prop2DAcoIsoDenQ_DEO2_FDTD -- nx,nz; $nx,$nz")
end
