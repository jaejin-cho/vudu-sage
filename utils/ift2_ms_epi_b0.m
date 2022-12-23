function [out] = ift2_ms_epi_b0( data, FTappa)

    [nx,ny,nc,ns] = size(data);

    out     =   zeross([nx,ny,nc,ns]);

    for xx = 1:nx
        for sh = 1:ns
            out(xx,:,:,sh) = ctranspose(sq(FTappa(:,:,xx,sh))) * reshape(data(xx,:,:,sh),[ny,nc]);
        end
    end

end
