#ifndef RTM_SEISMIC_H
#define RTM_SEISMIC_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <RTMBase.hpp>
#include <RTMGrid.hpp>
#include <RTMParam.hpp>
#include <RTMException.hpp>

/**
 * @brief Class RTMInstrumentRecord Inheritance of class Grid1D<GridData_type>
 */
template <typename GridData_type, typename DevPtr_type>
class RTMInstrumentRecord : public RTMVector<GridData_type, DevPtr_type>
{
protected:
    RTMGridCoordinate position; ///< 3D position in the grid
    size_t nsamples;///< number of samples / trace
public:
    /**
     * @brief Construct a new RTMInstrumentRecord object
     */
    RTMInstrumentRecord()
        : RTMVector<GridData_type, DevPtr_type>{}
    {
        nsamples = 0;
    }
    /**
     * @brief Construct a new RTMInstrumentRecord object
     */
    RTMInstrumentRecord(size_t x, size_t y, size_t z, size_t _nsamples)
        : RTMVector<GridData_type, DevPtr_type>{_nsamples}
    {
        position.setCoordinates(x, y, z);
        nsamples = _nsamples;
    }
    RTMInstrumentRecord(size_t x, size_t y, size_t z)
        : RTMVector<GridData_type, DevPtr_type>{}
    {
        position.setCoordinates(x, y, z);
    }
    /**
     * @brief Destroy the RTMInstrumentRecord object
     */
    ~RTMInstrumentRecord()
    {
    }
    /**
     * @brief Construct a new RTMInstrumentRecord object
     * @param v 
     */
    RTMInstrumentRecord(const RTMInstrumentRecord<GridData_type, DevPtr_type> &v)
        : RTMVector<GridData_type, DevPtr_type>{v}
    {
        // position.setCoordinates(v.getX(), v.getY(), v.getZ());
        position = v.position;
    }
    ///< Get methods
    size_t getX()
    {
        return (position.X);
    }
    size_t getY()
    {
        return (position.Y);
    }
    size_t getZ()
    {
        return (position.Z);
    }
    ///< Set methos
    void setX(size_t k)
    {
        position.X = k;
    }
    void setY(size_t k)
    {
        position.Y = k;
    }
    void setZ(size_t k)
    {
        position.Z = k;
    }
    /**
     * @brief Get the Position object
     * @return RTMGridCoordinate& 
     */
    RTMGridCoordinate &getPosition()
    {
        return position;
    }
};

/**
 * @brief Class RTMSeismicSource Inheritance of class RTMInstrumentRecord<GridData_type>
 */
template <typename GridData_type, typename DevPtr_type>
class RTMSeismicSource : public RTMInstrumentRecord<GridData_type, DevPtr_type>
{

protected:
    float dt;///<sampling interval in t
    float fpeak; ///< source peak frequency
public:
    /**
     * @brief Construct a new RTMSeismicSource object
     */
    RTMSeismicSource(size_t x, size_t y, size_t z, size_t nt)
        : RTMInstrumentRecord<GridData_type, DevPtr_type>{x, y, z}
    {
       RTMInstrumentRecord<GridData_type, DevPtr_type>::nsamples = nt;
    }
    /**
     * @brief Construct a new RTMSeismicSource object
     * @param x 
     * @param y 
     * @param z 
     * @param nt 
     * @param _dt 
     * @param _fpeak 
     */
    RTMSeismicSource(int x, int y, int z, int nt, float _dt, float _fpeak)
        : RTMInstrumentRecord<GridData_type, DevPtr_type>{x, y, z, nt}
    {
        dt = _dt;
        fpeak = _fpeak;
        RTMInstrumentRecord<GridData_type, DevPtr_type>::nsamples = nt;
        loadSamples(_dt, _fpeak);
    }
    /**
     * @brief Function loadSamples()
     */
    void loadSamples()
    {
        size_t it;
        RTMSeismicSource<GridData_type, DevPtr_type>::initGrid(RTMSeismicSource<GridData_type, DevPtr_type>::nsamples);
        for (it = 0; it < RTMSeismicSource<GridData_type, DevPtr_type>::size(); it++)
        {
            RTMSeismicSource<GridData_type, DevPtr_type>::set(static_cast<GridData_type>(ricker(it * dt - 1.0 / fpeak, fpeak)), it);
        }
        RTMSeismicSource<GridData_type, DevPtr_type>::updateHostPtr();
    }
    /**
     * @brief Function loadSamples()
     * @param _dt 
     * @param _fp 
     */
    void loadSamples(float _dt, float _fp){
        setDt(_dt);
        setFpeak(_fp);
        loadSamples();
    }
    /**
     * @brief Function unloadSamples()
     * 
     */
    void unloadSamples(){
        RTMInstrumentRecord<GridData_type, DevPtr_type>::destroyGrid();
    }
    /**
     * @brief Set the Dt object
     * @param _dt 
     */
    void setDt(float _dt){
        dt = _dt;
    }
    /**
     * @brief Set the Fpeak object
     * @param _fp 
     */
    void setFpeak(float _fp){
        fpeak=_fp;
    }
private:
    /**
     * @brief Function ricker()
     * @details Compute Ricker wavelet as a function of time.
     *          Notes: The amplitude of the Ricker wavelet at a frequency of 2.5*fpeak is 
     *          approximately 4 percent of that at the dominant frequency fpeak.
     *          The Ricker wavelet effectively begins at time t = -1.0/fpeak.  Therefore,
     *          for practical purposes, a causal wavelet may be obtained by a time delay
     *          of 1.0/fpeak.
     *          The Ricker wavelet has the shape of the second derivative of a Gaussian.
     * @author Dave Hale, Colorado School of Mines, 04/29/90
     * @param t time at which to evaluate Ricker wavelet
     * @param fpeak peak (dominant) frequency of wavelet
     * @return float 
     */
    float ricker(float t, float fpeak)
    {
        float x, xx;

        x = M_PI * fpeak * t;
        xx = x * x;
        return exp(-xx) * (1.0 - 2.0 * xx);
    }
};

/**
 * @brief Class RTMReceiverGrid Inheritance of class  RTMCube<GridData_type>
 * @tparam T 
 */
template <typename GridData_type, typename DevPtr_type>
class RTMReceiverGrid : public RTMCube<GridData_type, DevPtr_type>{

protected:
    size_t nt;///< number total of time steps
    size_t offsetX;
    size_t offsetY;
    size_t offsetZ;

    size_t distanceX;
    size_t distanceY;

    string fileName;
public:
    /**
     * @brief Construct a new RTMReceiverGrid object
     */
    RTMReceiverGrid(size_t nx, size_t ny, size_t _nt)
    :RTMCube<GridData_type, DevPtr_type>{nx, ny, _nt}
    {
        nt = _nt;
        offsetX = 0;
        offsetY = 0;
        offsetZ = 0;
        distanceX = 1;
        distanceY = 1;
    }
    
    ///< Get and set methods
    void setFileName(string &_str){
        fileName = _str;
    }
    string& getFileName(){
        return fileName;
    }
    void setDistanceX(size_t _k){
        distanceX = _k;
    }
    size_t getDistanceX(){
        return distanceX;
    }
    void setDistanceY(size_t _k){
        distanceY = _k;
    }
    size_t getDistanceY(){
        return distanceY;
    }
    void setOffsetX(size_t _k){
        offsetX = _k;
    }
    size_t getOffsetX(){
        return offsetX;
    }
    void setOffsetY(size_t _k){
        offsetY = _k;
    }
    size_t getOffsetY(){
        return offsetY;
    }
    void setOffsetZ(size_t _k){
        offsetZ = _k;
    }
    size_t getOffsetZ(){
        return offsetZ;
    }
};

/**
 * @brief Class RTMTaperFunction Inheritance of class Grid1D<GridData_type>
 * @tparam T 
 */
template <typename GridData_type, typename DevPtr_type>
class RTMTaperFunction : public RTMVector<GridData_type, DevPtr_type>
{

protected:
    size_t tplen; ///< taper length
    float tpfreq; // taper function frequency
public:
    RTMTaperFunction(size_t len, float freq)
        : RTMVector<GridData_type, DevPtr_type>{len}
    {
        tplen = len;
        tpfreq = freq;
        inittaper();
    }
    /**
     * @brief Function inittaper()
     */
    void inittaper()
    {
        size_t i;
        float dfrac;

        dfrac = sqrt(-log(tpfreq)) / (1. * tplen);
        /* taper vector is generated in decreasing order*/
        for (i = 0; i < tplen; i++)
        {
            RTMVector<GridData_type, DevPtr_type>::set(static_cast<GridData_type>(exp(-pow((dfrac * (i)), 2))), i);
        }
    }
    size_t tpLen()
    {
        return tplen;
    }
};

/**
 * @brief Class RTMShotDescriptor
 * @details Shot descriptor
 */
template <typename GridData_type, typename DevPtr_type>
class RTMShotDescriptor
{
protected:
    /* Each shot has an associated source record */
    RTMSeismicSource<GridData_type, DevPtr_type> *source = nullptr;

    /* Each shot has an associated list of receivers */
    RTMReceiverGrid<GridData_type, DevPtr_type> *recvGrid = nullptr;

    /* Shot migration image is stored here */
    RTMCube<GridData_type, DevPtr_type> *shotImage = nullptr;

    int NT;
public:
    /**
     * @brief Construct a new RTMShotDescriptor object
     */
    RTMShotDescriptor()
    {
    }
    /**
     * @brief Create a Shot Image object
     * @param nx 
     * @param ny 
     * @param nz 
     */
    void createShotImage(size_t nx, size_t ny, size_t nz)
    {
        shotImage = new RTMCube<GridData_type, DevPtr_type>(nx, ny, nz);
    }
    /**
     * @brief Destroy a Shot Image object
     * 
     */
    void destroyShotImage()
    {
        if (shotImage!=nullptr) delete shotImage;
    }
    ///< Get and set methods
    RTMCube<GridData_type, DevPtr_type>& getShotImage()
    {
        return *shotImage;
    }

    void setSource(RTMSeismicSource<GridData_type, DevPtr_type> *_source)
    {
        source = _source;
    }
    RTMSeismicSource<GridData_type, DevPtr_type> *getSource()
    {
        return source;
    }
    RTMReceiverGrid<GridData_type, DevPtr_type> *getReceiverGrid()
    {
        return recvGrid;
    }
    void setReceiverGrid(RTMReceiverGrid<GridData_type, DevPtr_type> &_recvGrid){
        recvGrid = _recvGrid;
    }

    int getNT()
    {
        return NT;
    }
    void setNT(int _nt)
    {
        NT = _nt;
    }

    /**
     * @brief Function loadReceiverGrid()
     * @details Load  receivers grid
     * @param rtmParam 
     */
    void loadReceiverGrid(RTMParam * rtmParam);
    /**
     * @brief Function unloadReceiverGrid()
     * @details Unload receivers grid
     */
    void unloadReceiverGrid(){
        if(recvGrid!=nullptr)
            delete recvGrid;
    }
    /**
     * @brief Function loadSourceSamples()
     * @details Load source samples vector
     */
    void loadSourceSamples(){
        if(source!=nullptr)
            source->loadSamples();
    }
    /**
     * @brief Function unloadSource()
     * @details Unload source samples vector
     */
    void unloadSource(){
        if(source!=nullptr){
            source->unloadSamples();
            delete source;
        }
    }
};

#endif