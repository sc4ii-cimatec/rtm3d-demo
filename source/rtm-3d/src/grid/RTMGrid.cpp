#include <iostream>
#include<vector> // for vector 
#include<algorithm> // for copy() and assign() 
#include<iterator> // for back_inserter 
#include <Misc.hpp>
#include <RTM.hpp>

using namespace std;

template <>
RTMData_t RTMCube<RTMData_t, RTMDevPtr_t>::getWithGlobalCoordinate(size_t gx, size_t gy, size_t gz)
{
    if( gx >= startX && gx < endX ){
        if( gy >= startY && gy < endY ){
            if( gz >= startZ && gz < endZ ){
                return get(gx-startX, gy-startY, gz-startZ);
            }
        }
    }
    return 0;
}
template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::setWithGlobalCoordinate(size_t gx, size_t gy, size_t gz, RTMData_t val)
{
    if( gx >= startX && gx < endX ){
        if( gy >= startY && gy < endY ){
            if( gz >= startZ && gz < endZ ){
                set(val, gx-startX, gy-startY, gz-startZ);
            }
        }
    }

}

template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::applyEnergy(size_t sx, size_t sy, size_t sz, RTMData_t _energy)
{
    RTMData_t val = getWithGlobalCoordinate(sx, sy, sz);
    setWithGlobalCoordinate(sx, sy, sz, (val + _energy));
}

/**
 * 
 * extends current grid borders
 * */
template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::extendBorders(uint32_t blength)
{

    if (blength > 0)
    {
        blen = blength;
    }
    this->setMaxMin(); // sets max and min values before extendind
    
    size_t ix, iy, iz, kx, ky, kz;
    size_t nx = getNX();
    size_t ny = getNY();
    size_t nz = getNZ();
    size_t nxe = nx + 2 * blen;
    size_t nye = ny + 2 * blen;
    size_t nze = nz + 2 * blen;
    HostBuffer_t<RTMData_t> tmpGrid;
    copy(getGridBuffer().begin(), getGridBuffer().end(), back_inserter(tmpGrid));
    
    
    size_t new_length = nxe*nye*nze;
    resize(nxe, nye, nze);
    for (ix = blen, kx = 0; ix < nxe - blen; ix++, kx++)
    {
        for (iy = blen, ky = 0; iy < nye - blen; iy++, ky++)
        {
            for (iz = blen, kz = 0; iz < nze - blen; iz++, kz++)
            {
                size_t offset_old;
                offset_old = kx*(ny*nz)+ (ky*nz) + kz;
                set(tmpGrid.at(offset_old), ix, iy, iz);
            }
        }
    }
    setNX(nxe);
    setNY(nye);
    setNZ(nze);
    // update device grid host pointer and length
    updateHostPtr();
}

template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::removeBorders(uint32_t blength)
{
    
    size_t ix, iy, iz, kx, ky, kz;
    size_t nx = getNX();
    size_t ny = getNY();
    size_t nz = getNZ();
    size_t nxr = getNX() - 2*blength;
    size_t nyr = getNY() - 2*blength;
    size_t nzr = getNZ() - 2*blength;
    HostBuffer_t<RTMData_t> tmpGrid;
    copy(getGridBuffer().begin(), getGridBuffer().end(), back_inserter(tmpGrid));
    
    size_t new_length = nxr*nyr*nzr;
    resize(nxr, nyr, nzr);
    for (ix = 0, kx = blength; ix < nxr; ix++, kx++)
    {
        for (iy = 0, ky = blength; iy < nyr; iy++, ky++)
        {
            for (iz = 0, kz = blength; iz < nzr; iz++, kz++)
            {
                size_t offset0;
                offset0 = (kx * (ny * nz)) + (ky * (nz)) + kz;
                set(tmpGrid.at(offset0), ix, iy, iz);
            }
        }
    }

    blen = blen - blength;
    if(blen < 0 ) blen=0;
    setNX(nxr);
    setNY(nyr);
    setNZ(nzr);
    setMaxMin(); // sets max and min values after removing

    // update device grid host pointer and length
    updateHostPtr();
}

template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::stack(RTMCube<RTMData_t, RTMDevPtr_t> &_rtmgrid)
{

    if (getNX() != _rtmgrid.getNX() ||
        getNY() != _rtmgrid.getNY() ||
        getNZ() != _rtmgrid.getNZ())
    {
        string msg = "[ RTMGrid::stack - Grids must be the same size! ]";
        GridException ex(msg);
        throw ex; 
    }
    int ix, iy, iz;
    for (ix = 0; ix < getNX(); ix++)
    {
        for (iy = 0; iy < getNY(); iy++)
        {
            for (iz = 0; iz < getNZ(); iz++)
            {
                RTMData_t val = get(ix, iy, iz) + _rtmgrid.get(ix, iy, iz);
                set(val, ix, iy, iz);
            }
        }
    }
}

template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::stackRegion(RTMCube<RTMData_t, RTMDevPtr_t> &regionGrid, int stX, int eX, int stY, int eY, int stZ, int eZ)
{

    assert((stX>=0));
    assert((stY>=0));
    assert((stZ>=0));
    assert((eX<=getNX()));
    assert((eY<=getNY()));
    assert((eZ<=getNZ()));
    assert((eX>stX));
    assert((eY>stY));
    assert((eZ>stZ));

    int lx = eX-stX;
    int ly = eY-stY;
    int lz = eZ-stZ;
    int ix, iy, iz;
    for (ix = stX; ix < eX; ix++)
    {
        for (iy = stY; iy < eY; iy++)
        {
            for (iz = stZ; iz < eZ; iz++)
            {
                RTMData_t val = get(ix, iy, iz) + regionGrid.get(ix-stX, iy-stY, iz-stZ);
                set(val, ix, iy, iz);
            }
        }
    }
}

template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::filter(RTMStencil<RTMData_t,RTMDevPtr_t, RTM_NDIM_3D> &_kernel)
{
    size_t ix=0, iy=0, iz=0;
    size_t st_order = _kernel.getOrder();
    size_t half_order = st_order / 2;
    size_t length = getNX()*getNY()*getNZ();
    RTMData_t *gauss = new RTMData_t[length];
    if (getNX()==1)
    { // 2D scenario
        #pragma omp parallel for private(iy, iz) collapse(2)
        for (iy = half_order; iy < (getNY() - half_order); iy++)
        {
            for (iz = half_order; iz < (getNZ() - half_order); iz++)
            {
                RTMData_t acmy = 0.0,acmz = 0.0,lapl = 0.0;
                size_t ky = 0, kz = 0;
                for (size_t io = 0; io <= st_order; io++)
                {
                    ky = iy - half_order + io;
                    kz = iz - half_order + io;

                    acmy += get(0, ky, iz) * _kernel.getStencilCoef(RTMDim::Ydim, io);
                    acmz += get(0, iy, kz) * _kernel.getStencilCoef(RTMDim::Zdim, io);
                }
                lapl = acmy + acmz;
                gauss[getOffset(0, iy, iz)]= lapl;
            }
        }
    }
    else
    { // 3D scenario
        #pragma omp parallel for private(ix, iy, iz) collapse(3)
        for (ix = half_order; ix < (getNX() - half_order); ix++)
        {
            for (iy = half_order; iy < (getNY() - half_order); iy++)
            {
                for (iz = half_order; iz < (getNZ() - half_order); iz++)
                {
                    RTMData_t acmx = 0.0,acmy = 0.0,acmz = 0.0,lapl = 0.0;
                    size_t kx = 0, ky = 0, kz = 0;
                    for (size_t io = 0; io <= st_order; io++)
                    {
                        kx = ix - half_order + io;
                        ky = iy - half_order + io;
                        kz = iz - half_order + io;

                        acmx += get(kx, iy, iz) * _kernel.getStencilCoef(RTMDim::Xdim, io);
                        acmy += get(ix, ky, iz) * _kernel.getStencilCoef(RTMDim::Ydim, io);
                        acmz += get(ix, iy, kz) * _kernel.getStencilCoef(RTMDim::Zdim, io);
                    }
                    lapl = acmx + acmy + acmz;
                    gauss[getOffset(ix, iy, iz)]= lapl;
                }
            }
        }
    }
    for (iy = 0; iy < (getNY()); iy++)
    {
        for (iz = 0; iz < (getNZ()); iz++)
        {
            if (getNX() == 1)
            {
                set(gauss[getOffset(0, iy, iz)], 0, iy, iz);
            }
            else
            {
                for (ix = 0; ix < (getNX()); ix++)
                {
                    set(gauss[getOffset(ix, iy, iz)], ix, iy, iz);
                }
            }
        }
    }
    delete gauss;
}

template <>
RTMPlane<RTMData_t,RTMDevPtr_t> *RTMCube<RTMData_t, RTMDevPtr_t>::get2DSection(size_t k, RTMDim secDim)
{
    
    string s("Trying to create 2D section out of grid limits.");
    RTMException ex(s);
    unsigned int n1 = 0, n2 = 0;
    size_t nx = getNX();
    size_t ny = getNY();
    size_t nz = getNZ();
    switch (secDim)
    {
    case RTMDim::Xdim:
        n1 = ny;
        n2 = nz;
        if(k<startX || k>=endX){
            throw ex; 
        }
        break;
    case RTMDim::Ydim:
        n1 = nx;
        n2 = nz;
        if(k<startY || k>=endY){
            throw ex; 
        }
        break;
    case RTMDim::Zdim:
        n1 = nx;
        n2 = ny;
        if(k<startZ || k>=endZ){
            throw ex; 
        }
        break;
    default:
        n1 = ny;
        n2 = nz;
        break;
    }

    RTMPlane<RTMData_t,RTMDevPtr_t> *secGrid = new RTMPlane<RTMData_t,RTMDevPtr_t>(n1, n2);

    for (int i1 = 0; i1 < n1; i1++)
    {
        for (int i2 = 0; i2 < n2; i2++)
        {

            RTMData_t val = 0.0;
            switch (secDim)
            {
            case RTMDim::Xdim:
                val = get(k, i1, i2);
                break;
            case RTMDim::Ydim:
                val = get(i1, k, i2);
                break;
            case RTMDim::Zdim:
                val = get(i1, i2, k);
                break;
            default:
                val = get(k, i1, i2);
                break;
            }
            secGrid->set(val,i1, i2);
        }
    }
    return secGrid;
}

template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::initBorders_ABC()
{

    unsigned int offset0, offset1, ix, iy, iz, kx, ky, kz;
    unsigned int lx, ly, lz, nx, ny, nz, blen;

    nx = getNX();
    ny = getNY();
    nz = getNZ();
    blen = getBorderLength();

    lx = nx - 2 * blen;
    ly = ny - 2 * blen;
    lz = nz - 2 * blen;

    for (ix = 0; ix < (nx); ix++)
    {
        for (iy = 0; iy < (ny); iy++)
        {
            for (iz = 0; iz < (nz); iz++)
            {
                float v = 0;
                if (ix <= blen)
                {
                    kx = blen;
                }
                else if (ix >= (lx + blen))
                {
                    kx = lx + blen - 1;
                }
                else
                {
                    kx = ix;
                }
                if (iy <= blen)
                {
                    ky = blen;
                }
                else if (iy >= (ly + blen))
                {
                    ky = ly + blen - 1;
                }
                else
                {
                    ky = iy;
                }
                if (iz <= blen)
                {
                    kz = blen;
                }
                else if (iz >= (lz + blen))
                {
                    kz = lz + blen - 1;
                }
                else
                {
                    kz = iz;
                }
                RTMData_t val = get(kx, ky, kz);
                set(val,ix, iy, iz);
            }
        }
    }
}

template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::initBorders_RBC()
{    
    unsigned int offset0, offset1, ix, iy, iz, kx, ky, kz, dist;
    unsigned int lx, ly, lz, nx, ny, nz, blen;
    RTMData_t val0;

    nx = getNX();
    ny = getNY();
    nz = getNZ();
    blen = getBorderLength();

    lx = nx - 2 * blen;
    ly = ny - 2 * blen;
    lz = nz - 2 * blen;

    // top and bottom cube
    #pragma omp parallel for private(ix, iy, iz, dist, val0) collapse(3)
    for (ix = blen; ix < (nx - blen); ix++)
    {
        for (iy = blen; iy < (ny - blen); iy++)
        {
            for (iz = 0; iz < blen; iz++)
            {
                
                // top
                dist = blen - iz;
                val0 = get(ix,iy,blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, dist, blen);
                set(val0,ix,iy,iz);

                // bottom
                val0 = get(ix,iy,(nz-blen-1));
                val0 = randbetween(val0, MAXVAL, MINVAL, dist, blen);
                set(val0,ix,iy,(nz-1)-iz);
            }
        }
    }

    // left and right cubes
    #pragma omp parallel for private(ix, iy, iz, dist, val0) collapse(3)
    for (ix = blen; ix < nx-blen; ix++)
    {
        for (iz = blen; iz < (nz - blen); iz++)
        {
            for (iy = 0; iy < blen; iy++)
            {
                // left
                val0 = get(ix,blen,iz);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen - iy, blen);
                set(val0,ix,iy,iz);

                // right
                val0 = get(ix,(ny-blen-1), iz);
                val0 = randbetween(val0, MAXVAL, MINVAL, iy, blen);
                set(val0,ix,(ny-blen)+iy,iz);
            }
        }
    }

    // front and back cubes
    #pragma omp parallel for private(ix, iy, iz, dist, val0) collapse(3)
    for (iy = blen; iy < ny-blen; iy++)
    {
        for (iz = blen; iz < (nz - blen); iz++)
        {
            for (ix = 0; ix < blen; ix++)
            {
                // front
                val0 = get(blen,iy,iz);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen - ix, blen);
                set(val0,ix,iy,iz);

                // back
                val0 = get((nx-blen-1),iy, iz);
                val0 = randbetween(val0, MAXVAL, MINVAL, ix, blen);
                set(val0, (nx-blen)+ix,iy,iz);
            }
        }
    }
    // now the corners...
    #pragma omp parallel for private(ix, iy, iz, dist, val0, kx, ky, kz) collapse(3)
    for (ix = 0; ix < blen; ix++){
        for (iz = 0; iz < blen; iz++){
            for (iy=0; iy<ny; iy++){
                ky = iy;
                if (iy<blen){
                    ky = blen;
                }else if (iy>=(ny-blen)){
                    ky = ny - blen -1;
                }
                
                // front-top
                val0 = get(blen,ky,blen);
                if (iy<blen){
                    val0 = randbetween(val0, MAXVAL, MINVAL, blen-iy, blen);
                }else if (iy>=(ny-blen)){
                    val0 = randbetween(val0, MAXVAL, MINVAL, iy-(ny-blen), blen);;
                }
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-ix, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-iz, blen);
                set(val0,ix,iy,iz);

                // back-top
                val0 = get(nx-blen-1,ky,blen);
                if (iy<blen){
                    val0 = randbetween(val0, MAXVAL, MINVAL, blen-iy, blen);
                }else if (iy>=(ny-blen)){
                    val0 = randbetween(val0, MAXVAL, MINVAL, iy-(ny-blen), blen);;
                }
                val0 = randbetween(val0, MAXVAL, MINVAL, ix, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-iz, blen);
                set(val0, (nx-blen)+ix,iy,iz);


                // front-bottom
                val0 = get(blen,ky,nz-blen-1);
                if (iy<blen){
                    val0 = randbetween(val0, MAXVAL, MINVAL, blen-iy, blen);
                }else if (iy>=(ny-blen)){
                    val0 = randbetween(val0, MAXVAL, MINVAL, iy-(ny-blen), blen);;
                }
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-ix, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, iz, blen);
                set(val0, ix,iy,(nz-blen)+iz);

                // back-bottom
                val0 = get(nx-blen-1,ky,nz-blen-1);
                if (iy<blen){
                    val0 = randbetween(val0, MAXVAL, MINVAL, blen-iy, blen);
                }else if (iy>=(ny-blen)){
                    val0 = randbetween(val0, MAXVAL, MINVAL, iy-(ny-blen), blen);;
                }
                val0 = randbetween(val0, MAXVAL, MINVAL, ix, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, iz, blen);
                set(val0, (nx-blen)+ix,iy,(nz-blen)+iz);
            }
        }
    }
    #pragma omp parallel for private(ix, iy, iz, dist, val0, kx, ky, kz) collapse(3)
    for (iy = 0; iy < blen; iy++){
        for (iz = 0; iz < blen; iz++){
            for (ix=blen; ix<nx-blen; ix++){
                // left-top
                val0 = get(ix, blen, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-iy, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-iz, blen);
                set(val0,ix,iy,iz);

                // right-top
                val0 = get(ix, ny-blen-1, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, iy, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-iz, blen);
                set(val0,ix,(ny-blen)+iy,iz);

                // left-bottom
                val0 = get(ix, blen, nz-blen-1);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-iy, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, iz, blen);
                set(val0, ix,iy,(nz-blen)+iz);

                // right-bottom
                val0 = get(ix, ny-blen-1, nz-blen-1);
                val0 = randbetween(val0, MAXVAL, MINVAL, iy, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, iz, blen);
                set(val0, ix,(ny-blen)+iy,(nz-blen)+iz);
            }
        }
    }
    #pragma omp parallel for private(ix, iy, iz, dist, val0, kx, ky, kz) collapse(3)
    for (ix = 0; ix < blen; ix++){
        for (iy = 0; iy < blen; iy++){
            for (iz=blen; iz<nz-blen; iz++){
                // front-left
                val0 = get(blen, blen, iz);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-ix, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-iy, blen);
                set(val0,ix,iy,iz);

                // front-right
                val0 = get(blen, ny-blen-1, iz);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-ix, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, iy, blen);
                set(val0,ix,(ny-blen)+iy,iz);

                // back-left
                val0 = get(nx-blen-1, blen, iz);
                val0 = randbetween(val0, MAXVAL, MINVAL, ix, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, blen-iy, blen);
                set(val0, (nx-blen)+ix,iy,iz);

                // back-right
                val0 = get(nx-blen-1, ny-blen-1, iz);
                val0 = randbetween(val0, MAXVAL, MINVAL, ix, blen);
                val0 = randbetween(val0, MAXVAL, MINVAL, iy, blen);
                set(val0, (nx-blen)+ix,(ny-blen)+iy,iz);
            }
        }
    } 
}

template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::initBorders_HBC()
{
    unsigned int offset0, offset1, ix, iy, iz, kx, ky, kz;
    unsigned int lx, ly, lz, nx, ny, nz, blen;
    RTMData_t val0, val1, MAXVAL, MINVAL;

    nx = getNX();
    ny = getNY();
    nz = getNZ();
    blen = getBorderLength();

    MAXVAL = MAXVAL;
    MINVAL = MINVAL;

    lx = nx - 2 * blen;
    ly = ny - 2 * blen;
    lz = nz - 2 * blen;

    initBorders_RBC();

    /* Upper cube with ABC (Dussaud 2008) */
    for (ix = 0; ix < (nx); ix++)
    {
        for (iy = 0; iy < (ny); iy++)
        {
            for (iz = 0; iz < blen; iz++)
            {
                if (ix <= blen)
                {
                    kx = blen;
                }
                else if (ix >= (lx + blen))
                {
                    kx = lx + blen - 1;
                }
                else
                {
                    kx = ix;
                }
                if (iy <= blen)
                {
                    ky = blen;
                }
                else if (iy >= (ly + blen))
                {
                    ky = ly + blen - 1;
                }
                else
                {
                    ky = iy;
                }
                kz = (blen);

                set(get(kx, ky, kz), ix, iy, iz);
            }
        }
    }
}

template <>
void RTMCube<RTMData_t, RTMDevPtr_t>::initBorders(const RTMBoundaryCondition bc)
{

    switch (bc)
    {
    case RTMBoundaryCondition::ABC:
        initBorders_ABC();
        break;
    case RTMBoundaryCondition::HBC:
        initBorders_HBC();
        break;
    case RTMBoundaryCondition::RBC:
        initBorders_RBC();
        break;
    default:
        break;
    }
}

template<>
RTMCube<RTMData_t, RTMDevPtr_t> * RTMCube<RTMData_t, RTMDevPtr_t>::getSubGrid(size_t stX, size_t eX, size_t stY, size_t eY, size_t stZ, size_t eZ){


    assert((stX>=0));
    assert((stY>=0));
    assert((stZ>=0));
    assert((eX<=getNX()));
    assert((eY<=getNY()));
    assert((eZ<=getNZ()));
    assert((eX>stX));
    assert((eY>stY));
    assert((eZ>stZ));

    int lx = eX-stX;
    int ly = eY-stY;
    int lz = eZ-stZ;
    RTMCube<RTMData_t, RTMDevPtr_t> *ngrid = new RTMCube<RTMData_t, RTMDevPtr_t>(lx, ly,lz);
    ngrid->setStartX(stX);
    ngrid->setStartY(stY);
    ngrid->setStartZ(stZ);
    ngrid->setEndX(eX);
    ngrid->setEndY(eY);
    ngrid->setEndZ(eZ);
    int ix, iy, iz;
    for (ix=stX; ix<eX; ix++){
        for (iy=stY; iy<eY; iy++){
            for (iz=stZ; iz<eZ; iz++){
                RTMData_t val0 = get(ix,iy,iz);
                ngrid->set(val0, ix-stX,iy-stY,iz-stZ);
            }
        }
    }
    return ngrid;
}

template<>
RTMVelocityModel<RTMData_t, RTMDevPtr_t> * RTMVelocityModel<RTMData_t, RTMDevPtr_t>::createDirectWaveFilterGrid(int depthZ)
{
    size_t ix, iy, iz;
    size_t nx = getNX();
    size_t ny = getNY();
    size_t nz = getNZ();
    size_t blen = getBorderLength();
    RTMVelocityModel<RTMData_t, RTMDevPtr_t> * fGrid = new RTMVelocityModel<RTMData_t, RTMDevPtr_t>(nx, ny, nz);
    fGrid->setBorderLength(blen);

    for (ix = 0; ix < (nx); ix++)
    {
        for (iy = 0; iy < (ny); iy++)
        {
            for (iz = 0; iz < nz; iz++)
            {
                fGrid->set(get(ix,iy,depthZ), ix,iy,iz);
            }
        }
    }
    return fGrid;
}