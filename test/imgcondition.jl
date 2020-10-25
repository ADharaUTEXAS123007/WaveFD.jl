using LinearAlgebra, Test, WaveFD, PyPlot

dz,dy,dx = 10.0,10.0,10.0
nz,ny,nx = 101,101,101
nbz,nby,nbx = 32,32,32
nthreads = Sys.CPU_THREADS
n1 = 6
n2 = 12
z0 = div(nz,2)+1
x0 = div(nx,2)+1

@testset "ImgCondition" begin
    @testset "Imaging Condition Iso 2D tests" begin
        pH = WaveFD.Prop2DAcoIsoDenQ_DEO2_FDTD(freesurface=false, nz=nz, nx=nx, nbz=nbz, nbx=nbx, dz=10.0, dx=10.0, dt=0.001, nthreads=nthreads)
        bH,vH = WaveFD.B(pH), WaveFD.V(pH)
        bH .= 1
        vH .= 1
        poldH = WaveFD.POld(pH)
        poldH .= 0
        poldH[z0-n1:z0+n1,:] .= +1
        wH = zeros(Float32,nz,nx)
        wH[z0-n2:z0+n2,:] .= +1

        pV = WaveFD.Prop2DAcoIsoDenQ_DEO2_FDTD(freesurface=false, nz=nz, nx=nx, nbz=nbz, nbx=nbx, dz=10.0, dx=10.0, dt=0.001, nthreads=nthreads)
        bV,vV = WaveFD.B(pV), WaveFD.V(pV)
        bV .= 1
        vV .= 1
        poldV = WaveFD.POld(pV)
        poldV .= 0
        poldV[:,x0-n1:x0+n1] .= +1
        wV = zeros(Float32,nz,nx)
        wV[:,x0-n2:x0+n2] .= +1

        dstdH = zeros(Float32,nz,nx)
        dfwiH = zeros(Float32,nz,nx)
        drtmH = zeros(Float32,nz,nx)
        dstdV = zeros(Float32,nz,nx)
        dfwiV = zeros(Float32,nz,nx)
        drtmV = zeros(Float32,nz,nx)

        WaveFD.adjointBornAccumulation!(pH,WaveFD.ImagingConditionStandard(),dstdH,wH)
        WaveFD.adjointBornAccumulation!(pH,WaveFD.ImagingConditionWaveFieldSeparationFWI(),dfwiH,wH)
        WaveFD.adjointBornAccumulation!(pH,WaveFD.ImagingConditionWaveFieldSeparationRTM(),drtmH,wH)

        WaveFD.adjointBornAccumulation!(pV,WaveFD.ImagingConditionStandard(),dstdV,wV)
        WaveFD.adjointBornAccumulation!(pV,WaveFD.ImagingConditionWaveFieldSeparationFWI(),dfwiV,wV)
        WaveFD.adjointBornAccumulation!(pV,WaveFD.ImagingConditionWaveFieldSeparationRTM(),drtmV,wV)

        println("")
        @show extrema(dstdH)
        @show extrema(dfwiH)
        @show extrema(drtmH)

        println("")
        @show extrema(dstdV)
        @show extrema(dfwiV)
        @show extrema(drtmV)

        figure(figsize=(8,6.5)); clf()
        subplot(2,2,1); imshow(poldV, cmap="seismic", aspect="auto"); clim(-2,+2); colorbar(); title("p Vertical")
        subplot(2,2,2); imshow(wV, cmap="seismic", aspect="auto"); clim(-2,+2); colorbar(); title("wavefield Vertical")
        subplot(2,2,3); imshow(poldH, cmap="seismic", aspect="auto"); clim(-2,+2); colorbar(); title("p Horizontal")
        subplot(2,2,4); imshow(wH, cmap="seismic", aspect="auto"); clim(-2,+2); colorbar(); title("wavefield Horizontal")
        tight_layout()
        display(gcf())
        savefig("imgcondition.input.png")

        figure(figsize=(8,6.5)); clf()
        subplot(2,2,1); imshow(dstdH, cmap="seismic", aspect="auto"); clim(-2,+2); colorbar(); title("Standard IC Horizontal")
        subplot(2,2,2); imshow(dfwiH, cmap="seismic", aspect="auto"); clim(-2,+2); colorbar(); title("FWI IC Horizontal")
        subplot(2,2,3); imshow(drtmH, cmap="seismic", aspect="auto"); clim(-2,+2); colorbar(); title("RTM IC Horizontal")
        tight_layout()
        display(gcf())
        savefig("imgcondition.output.horizontal.png")

        figure(figsize=(8,6.5)); clf()
        subplot(2,2,1); imshow(dstdV, cmap="seismic", aspect="auto"); clim(-2,+2); colorbar(); title("Standard IC Vertical")
        subplot(2,2,2); imshow(dfwiV, cmap="seismic", aspect="auto"); clim(-2,+2); colorbar(); title("FWI IC Vertical")
        subplot(2,2,3); imshow(drtmV, cmap="seismic", aspect="auto"); clim(-2,+2); colorbar(); title("RTM IC Vertical")
        tight_layout()
        display(gcf())
        savefig("imgcondition.output.vertical.png")
    end
end

nothing