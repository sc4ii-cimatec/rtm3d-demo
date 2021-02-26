#ifndef RTMGRID_H
#define RTMGRID_H
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <Misc.hpp>
#include <Grid.hpp>
#include <RTMBase.hpp>

using namespace std;
using namespace grid;

// RTMGrid definitions
#define DEFAULT_BORDER_SIZE         40
#define RTMGRID_MAX_DIM             4


/**
 * @brief Enum Class RTMBoundaryCondition
 */
enum class RTMBoundaryCondition
{
    NONE,
    ABC,
    HBC,
    RBC
};
/**
 * @brief Inline function getRTMBoundaryCondition()
 * @param str 
 * @return RTMBoundaryCondition 
 */
static inline RTMBoundaryCondition getRTMBoundaryCondition(string str)
{
    if (str.compare("hbc") == 0)
    {
        return RTMBoundaryCondition::HBC;
    }
    else if (str.compare("rbc") == 0)
    {
        return RTMBoundaryCondition::RBC;
    }
    else if (str.compare("abc") == 0)
    {
        return RTMBoundaryCondition::ABC;
    }
    else
    {
        return RTMBoundaryCondition::HBC;
    }
}
/**
 * @brief Enum class RTMDim
 */
enum class RTMDim
{
    Xdim = 0,
    Ydim = 1,
    Zdim = 2
};

/**
 * @brief Class RTMGridCoordinate
 */
class RTMGridCoordinate
{
public:
    int X, Y, Z;
    /**
     * @brief Construct a new RTMGridCoordinate object
     */
    RTMGridCoordinate()
    {
        X = 0;
        Y = 0;
        Z = 0;
    }
    /**
     * @brief Destroy the RTMGridCoordinate object
     */
    ~RTMGridCoordinate()
    {
    }
    /**
     * @brief Copy Construct a new RTMGridCoordinate object
     * @param v reference to RTMGridCoordinate
     */
    RTMGridCoordinate(const RTMGridCoordinate &v)
    {
        X = v.X;
        Y = v.Y;
        Z = v.Z;
    }
    /**
     * @brief Assignment operator
     * @param v Reference to RTMGridCoordinate
     * @return RTMGridCoordinate& 
     */
    RTMGridCoordinate &operator=(const RTMGridCoordinate &v)
    {
        X = v.X;
        Y = v.Y;
        Z = v.Z;
        return *this;
    }
    /**
     * @brief Construct a new RTMGridCoordinate object
     * @param x Integer value
     * @param y Integer value
     * @param z Integer value
     */
    RTMGridCoordinate(int x, int y, int z)
    {
        X = x;
        Y = y;
        Z = z;
    }
    /**
     * @brief Set the Coordinates object
     * @param x Integer value
     * @param y Integer value
     * @param z Integer value
     */
    void setCoordinates(int x, int y, int z)
    {
        X = x;
        Y = y;
        Z = z;
    }
    /**
     * @brief Get X
     * @return int 
     */
    int getX()
    {
        return X;
    }
    /**
     * @brief Get Y
     * @return int 
     */
    int getY()
    {
        return Y;
    }
    /**
     * @brief Get Z
     * @return int 
     */
    int getZ()
    {
        return Z;
    }
    /**
     * @brief Function toString()
     * @return string 
     */
    string toString()
    {
        std::ostringstream oss;
        oss << "(" << X << "," << Y << "," << Z << ")";
        return oss.str();
    }
};

template <typename GridData_type, typename DevPtr_type, unsigned int _nDim=1>
class RTMGrid : public Grid<GridData_type>, public DefaultDeviceGrid<GridData_type, DevPtr_type>
{
protected:
    uint32_t    nDim;
    size_t      dimSizes[RTMGRID_MAX_DIM];
public:
    RTMGrid()
    :Grid<GridData_type>{},DefaultDeviceGrid<GridData_type, DevPtr_type>{}
    {
        nDim=_nDim; 
        for(int i=0; i<RTMGRID_MAX_DIM; i++)dimSizes[i] = 1;
    }
    RTMGrid(size_t nt, size_t nx=1, size_t ny=1, size_t nz=1)
    :Grid<GridData_type>{nt*nx*ny*nz},DefaultDeviceGrid<GridData_type, DevPtr_type>{}
    {
        for(int i=0; i<RTMGRID_MAX_DIM; i++)dimSizes[i] = 1;
        nDim=_nDim;
        if(nDim>=1) dimSizes[0] = nt;
        if(nDim>=2) dimSizes[1] = nx;
        if(nDim>=3) dimSizes[2] = ny;
        if(nDim>=4) dimSizes[3] = nz;
        size_t len = nt*nx*ny*nz;
        DefaultDeviceGrid<GridData_type, DevPtr_type>::HOST_PTR = Grid<GridData_type>::data();
        DefaultDeviceGrid<GridData_type, DevPtr_type>::BUFFER_LENGTH = len*sizeof(GridData_type);
    }
    void reset(){
        Grid<GridData_type>::fill(0.0);
        DefaultDeviceGrid<GridData_type, DevPtr_type>::devMemSet(0.0);
    }

    void resize(size_t k0, size_t k1=0, size_t k2=0, size_t k3=0){

        size_t newsize = 1;
        if (k0!=0) {newsize*=k0; dimSizes[0] = k0;}
        if (k1!=0) {newsize*=k1; dimSizes[1] = k1;}
        if (k2!=0) {newsize*=k2; dimSizes[2] = k2;};
        if (k3!=0) {newsize*=k3; dimSizes[3] = k3;};

        if (k0!=0){
            Grid<GridData_type>::resizeGrid(newsize);
        }
    }


    size_t getOffset(size_t k0, size_t k1=0, size_t k2=0, size_t k3=0){
        size_t offset = (k0*dimSizes[3]*dimSizes[2]*dimSizes[1]) + 
        (k1*dimSizes[3]*dimSizes[2]) + (k2*dimSizes[3]) + k3;
        return offset;
    }

    GridData_type get(size_t k0, size_t k1=0, size_t k2=0, size_t k3=0)
    {
        size_t offset = getOffset(k0,k1,k2,k3);
        return Grid<GridData_type>::getByOffset(offset);
    }
    void set(GridData_type val, size_t k0, size_t k1=0, size_t k2=0, size_t k3=0)
    {
        size_t offset = getOffset(k0,k1,k2,k3);
        Grid<GridData_type>::setByOffset(offset, val);
    }
    void updateHostPtr()
    {   
        DefaultDeviceGrid<GridData_type, DevPtr_type>::HOST_PTR = Grid<GridData_type>::data();
        DefaultDeviceGrid<GridData_type, DevPtr_type>::BUFFER_LENGTH = Grid<GridData_type>::gridSize*sizeof(GridData_type); 
    }
    size_t getNX(){
        if(nDim >= 4){
            return dimSizes[1];
        }else{
            return dimSizes[0];
        }
    }
    size_t getNY(){
        if(nDim >= 4){
            return dimSizes[2];
        }else{
            return dimSizes[1];
        }
    }
    size_t getNZ(){
        if(nDim >= 4){
            return dimSizes[3];
        }else{
            return dimSizes[2];
        }
    }
    size_t getNT(){
        return dimSizes[0];
    }
protected:
    void setNX(size_t _k){
        if(nDim >= 4){
            dimSizes[1] = _k;
        }else{
            dimSizes[0] = _k;
        }
    }
    void setNY(size_t _k){
        if(nDim >= 4){
            dimSizes[2] = _k;
        }else{
            dimSizes[1] = _k;
        }
    }
    void setNZ(size_t _k){
        if(nDim >= 4){
            dimSizes[3] = _k;
        }else{
            dimSizes[2] = _k;
        }
    }
    void setNT(size_t _k){
        dimSizes[0] = _k;
    }
};

template <typename GridData_type, typename DevPtr_type>
class RTMVector : public RTMGrid<GridData_type, DevPtr_type, 1U>
{
public:
    RTMVector(size_t nx) 
    : RTMGrid<GridData_type, DevPtr_type, 1U>{nx}
    {}
    RTMVector() 
    : RTMGrid<GridData_type, DevPtr_type, 1U>{}
    {}
};

template <typename GridData_type, typename DevPtr_type>
class RTMPlane : public RTMGrid<GridData_type, DevPtr_type, 2U>
{
public:
     RTMPlane(size_t nx, size_t ny) 
    : RTMGrid<GridData_type, DevPtr_type, 2U>{nx, ny}
    {}
};

/**
 * @brief Class RTMFilterKernel Inheritance of class Grid2D<GridData_type>
 * @tparam GridData_type 
 */
template <typename GridData_type, typename DevPtr_type>
class RTMFilterKernel : public RTMVector<GridData_type, DevPtr_type>
{
protected:
    uint32_t order; ///< Filter order
public:
    RTMFilterKernel(uint32_t _order)
        : RTMVector<GridData_type, DevPtr_type>{_order + 1}
    {
        order = _order;
        setLaplacianKernel(); // Default Mode: Laplacian
    }

    RTMFilterKernel(const RTMFilterKernel<GridData_type, DevPtr_type> &_k)
        : RTMVector<GridData_type, DevPtr_type>{(_k.getOrder() + 1)}
    {
        order = _k.getOrder();
        copyData(_k);
    }
    /**
     * @brief Construct a new RTMFilterKernel object
     * @param _ndims Integer value
     * @param _order Integer value
     */
    RTMFilterKernel & operator=(const RTMFilterKernel<GridData_type, DevPtr_type> &_k) 
    {
        order = _k.getOrder();
        RTMVector<GridData_type, DevPtr_type>::initGrid(order+1);
        copyData(_k);
        return this;
    }

    void setLaplacianKernel(){
        laplacian_coefs(RTMVector<GridData_type, DevPtr_type>::data(), order);
    }

    /**
     * @brief Get the Order object
     * @return int 
     */
    int getOrder()
    {
        return order;
    }
    /**
     * @brief Get the Coef object
     * @param k 
     * @return float 
     */
    float getCoef(size_t k)
    {
        return RTMVector<GridData_type, DevPtr_type>::get(k);
    }

};
/**
 * @brief Class RTMStencil Inheritance of class RTMFilterKernel<GridData_type> 
 * @tparam GridData_type 
 */
template <typename GridData_type, typename DevPtr_type, uint32_t _nDims=3>
class RTMStencil : public RTMPlane<GridData_type, DevPtr_type>
{
protected:
    RTMFilterKernel<GridData_type, DevPtr_type> *kernel;
    uint32_t stencilOrder;
public:
    /**
     * @brief Construct a new RTMStencil object
     * @param _op 
     */
    RTMStencil(const RTMStencil<GridData_type, DevPtr_type> &_op)
        : RTMPlane<GridData_type, DevPtr_type>{_op}
    {
        assert(_nDims==_op.getNDims());
        stencilOrder=_op.getOrder();
        kernel = new RTMFilterKernel<GridData_type, DevPtr_type>(stencilOrder);
        kernel->copyData(_op.getKernel());
    }
    RTMStencil & operator=(const RTMStencil<GridData_type, DevPtr_type> &_op){
        assert(_nDims==_op.getNDims());
        stencilOrder=_op.getOrder();
        kernel = new RTMFilterKernel<GridData_type, DevPtr_type>(stencilOrder);
        kernel->copyData(_op.getKernel());
        RTMPlane<GridData_type, DevPtr_type>::copyData(_op);
    }

    RTMStencil(uint32_t _order)
        : RTMPlane<GridData_type, DevPtr_type>{_nDims, (_order+1)}
    {
        stencilOrder=_order;
        kernel = new RTMFilterKernel<GridData_type, DevPtr_type>(stencilOrder);
        kernel->setLaplacianKernel();
    }

    RTMStencil(uint32_t _order, GridData_type* derivatives)
        : RTMPlane<GridData_type, DevPtr_type>{_nDims, (_order+1)}
    {
        stencilOrder=_order;
        kernel = new RTMFilterKernel<GridData_type, DevPtr_type>(_order);
        kernel->setLaplacianKernel();
        setStencilCoefficients(derivatives);
    }
    GridData_type getStencilCoef(RTMDim dim, uint32_t k)
    {
        return RTMPlane<GridData_type, DevPtr_type>::get(uint32_t(dim), k);
    }
    
    /**
    * Computes and stores inverted squared spatial derivatives for laplacian calculation.
    * 
    * dx[k] = [(1.0/dx)^2]*coefs[k]
    * dy[k] = [(1.0/dy)^2]*coefs[k]
    * dz[k] = [(1.0/dz)^2]*coefs[k]
    * 
    * Derivatives must be passed inside a vector ndims long.
    * 
    * CONVENTION:  deriv[0] = X; deriv[1] = Y; deriv[2] = Z
    * 
    */
    void setStencilCoefficients(GridData_type* derivatives)
    {

        for (int d = 0; d < _nDims; d++)
        {
            float dev = derivatives[d];
            float devinv = (1.0 / dev) * (1.0 / dev);
            for (int io = 0; io < (stencilOrder + 1); io++)
            {
                float val = devinv * kernel->getCoef(io);
                RTMPlane<GridData_type, DevPtr_type>::set(val, d,io);
            }
        }
    }

    GridData_type* getStencilCoefArray(RTMDim dim)
    {
        GridData_type * arr = new GridData_type[stencilOrder + 1];
        int k;
        for (k=0; k<(stencilOrder + 1); k++){
            arr[k] = getStencilCoef(dim, k);
        }
        return arr;
    }
    RTMVector<GridData_type, DevPtr_type> & getStencilCoefVector(RTMDim dim){
        RTMVector<GridData_type, DevPtr_type> * arr = new RTMVector<GridData_type, DevPtr_type>(stencilOrder + 1);
        int k;
        for (k=0; k<(stencilOrder + 1); k++){
            arr->set(getStencilCoef(dim, k), k);
        }
        return *arr;
    }

    RTMFilterKernel<GridData_type, DevPtr_type> & getKernel(){return *kernel;}
    uint32_t getOrder(){return stencilOrder;}
    uint32_t getNDims(){return _nDims;}
};


/**
 * @brief Class RTMGrid Inheritance of class Grid3D<GridData_type>
 * @tparam GridData_type 
*/
template <typename GridData_type, typename DevPtr_type>
class RTMCube : public RTMGrid<GridData_type, DevPtr_type, 3U>
{
protected:
    string vfilename;
    uint32_t blen;

    ///< Processing coordinates
    size_t startX;
    size_t startY;
    size_t startZ;

    size_t endX;
    size_t endY;
    size_t endZ;

public:
    /**
     * @brief Construct a new RTMGrid object
     */
    RTMCube(size_t lx, size_t ly, size_t lz)
        : RTMGrid<GridData_type, DevPtr_type, 3U>{lx, ly, lz}
    {
        blen = DEFAULT_BORDER_SIZE;
        startX = 0;
        startY = 0;
        startZ = 0;
        endX = lx;
        endY = ly;
        endZ = lz;
    }
    /**
     * @brief Construct a new RTMGrid object
     */
    RTMCube(size_t lx, size_t ly, size_t lz, string fname)
        : RTMGrid<GridData_type, DevPtr_type, 3U>{lx, ly, lz}
    {
        blen = DEFAULT_BORDER_SIZE;
        vfilename = fname;
    }
    void setBorderLength(uint32_t b)
    {
        blen = b;
    }
    uint32_t getBorderLength() const
    {
        return blen;
    }
    /* Extends grid's borders in all directions */
    void extendBorders(uint32_t blength);
    /* Removes grid's borders in all directions */
    void removeBorders(uint32_t blength);
    /* Fills borders values according to chosen Boundary Condition method*/
    void initBorders(const RTMBoundaryCondition bc);
    /* Stacks another grid image onto current grid (+=)*/
    void stack(RTMCube<GridData_type,DevPtr_type> &_rtmgrid);
    /* Stacks another grid's region image onto the same region in current grid (+=)*/
    void stackRegion(RTMCube<GridData_type,DevPtr_type> &regionGrid, int stX, int eX, int stY, int eY, int stZ, int eZ);
    /* Filters current grid using kernel param */
    void filter(RTMStencil<RTMData_t,RTMDevPtr_t, RTM_NDIM_3D> &_kernel);
    /* Extract a 2D section at a given "pos" position in "secDim" dimension */
    RTMPlane<GridData_type,DevPtr_type> *get2DSection(size_t pos, RTMDim secDim);

    ///< Set methods
    void setStartX(size_t sX){startX = sX;};
    void setStartY(size_t sY){startY = sY;};
    void setStartZ(size_t sZ){startZ = sZ;};
    void setEndX(size_t eX){endX = eX;};
    void setEndY(size_t eY){endY = eY;};
    void setEndZ(size_t eZ){endZ = eZ;};
    ///< Get methods
    size_t getStartX(){return startX;};
    size_t getStartY(){return startY;};
    size_t getStartZ(){return startZ;};
    size_t getEndX(){return endX;};
    size_t getEndY(){return endY;};
    size_t getEndZ(){return endZ;};

    /**
     * These get and set functions receive xyz coordinates relative to the
     * global grid values, not local ones (Used in taper functions)
     * */
    GridData_type getWithGlobalCoordinate(size_t gx, size_t gy, size_t gz);
    void setWithGlobalCoordinate(size_t gx, size_t gy, size_t gz, GridData_type t);

    void applyEnergy(size_t x, size_t y, size_t z, GridData_type _energy);
    RTMCube<GridData_type, DevPtr_type> * getSubGrid(size_t startx, size_t endx, 
    size_t starty, size_t endy, size_t startz, size_t endz);

private:
    void initBorders_HBC();
    void initBorders_RBC();
    void initBorders_ABC();
};
/**
 * @brief Class RTMVelocityModel Inheritance of class  RTMCube<GridData_type>
 * @tparam T 
 */
template <typename GridData_type, typename DevPtr_type>
class RTMVelocityModel : public RTMCube<GridData_type, DevPtr_type>
{
protected:
    string modelname; ///< modeling name
public:
    /**
     * @brief Construct a new RTMVelocityModel object
     * @param lx Integer value
     * @param ly Integer value
     * @param lz Integer value
     */
    RTMVelocityModel(size_t lx, size_t ly, size_t lz)
        : RTMCube<GridData_type, DevPtr_type>(lx, ly, lz)
    {
    }
    /**
     * @brief Construct a new RTMVelocityModel object
     */
    RTMVelocityModel(size_t lx, size_t ly, size_t lz, string mname)
        : RTMCube<GridData_type, DevPtr_type>(lx, ly, lz)
    {
        modelname = mname;
    }
    /**
     * @brief Construct a new RTMVelocityModel object
     * @param mname: model name
     * @param filename : model file
     */
    RTMVelocityModel(size_t lx, size_t ly, size_t lz, string mname, string filename)
        : RTMCube<GridData_type, DevPtr_type>(lx, ly, lz)
    {
        RTMCube<GridData_type, DevPtr_type>::vfilename = filename;
        modelname = mname;
        RTMCube<GridData_type, DevPtr_type>::loadFromFile(filename);
        RTMCube<GridData_type, DevPtr_type>::setMaxMin();
    }

    RTMVelocityModel<GridData_type, DevPtr_type> * createDirectWaveFilterGrid(int depthZ);

};

/**
 * @brief Class RTMGridCollection Inheritance of class Grid4D<GridData_type>
 * @tparam T 
 */
template <typename GridData_type, typename DevPtr_type>
class RTMGridCollection : public RTMGrid<GridData_type, DevPtr_type, 4U>
{
/**
 * A RTMGrid Collection can be understood as a 
 * collection of 3D grids over time. For example,
 * to store a history of upper-border cubes (UPB) on HBC
 * migrations.
 */
protected:
    string vfilename;
    uint32_t blen;///< border length

    ///< Processing coordinates
    size_t startT;
    size_t startX;
    size_t startY;
    size_t startZ;

    size_t endT;
    size_t endX;
    size_t endY;
    size_t endZ;

public:
    /**
     * @brief Construct a new RTMGridCollection object
     */
    RTMGridCollection(size_t lt, size_t lx, size_t ly, size_t lz)
        : RTMGrid<GridData_type, DevPtr_type, 4U>{lt, lx, ly, lz}
    {
        blen = DEFAULT_BORDER_SIZE;
        startT = 0;
        startX = 0;
        startY = 0;
        startZ = 0;
        endT = lt;
        endX = lx;
        endY = ly;
        endZ = lz;
    }
    /**
     * @brief Construct a new RTMGridCollection object
     */
    RTMGridCollection(size_t lt, size_t lx, size_t ly, size_t lz, string fname)
        : RTMGrid<GridData_type, DevPtr_type, 4U>{lt, lx, ly, lz}
    {
        blen = DEFAULT_BORDER_SIZE;
        vfilename = fname;
        startT = 0;
        startX = 0;
        startY = 0;
        startZ = 0;
        endT = lt;
        endX = lx;
        endY = ly;
        endZ = lz;
    }
    /**
     * @brief Set the Border Length
     * @param b 
     */
    void setBorderLength(uint32_t b)
    {
        blen = b;
    }
    /**
     * @brief Get the Border Length
     */
    uint32_t getBorderLength()
    {
        return blen;
    }

    void setStartT(size_t sT){startT = sT;};
    void setStartX(size_t sX){startX = sX;};
    void setStartY(size_t sY){startY = sY;};
    void setStartZ(size_t sZ){startZ = sZ;};
    void setEndT(size_t eT){endT = eT;};
    void setEndX(size_t eX){endX = eX;};
    void setEndY(size_t eY){endY = eY;};
    void setEndZ(size_t eZ){endZ = eZ;};
    ///< Get methods
    size_t getStartT(){return startT;};
    size_t getStartX(){return startX;};
    size_t getStartY(){return startY;};
    size_t getStartZ(){return startZ;};
    size_t getEndT(){return endT;};
    size_t getEndX(){return endX;};
    size_t getEndY(){return endY;};
    size_t getEndZ(){return endZ;};
};

#endif
