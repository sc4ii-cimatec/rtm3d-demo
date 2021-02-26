#ifndef GRID_H
#define GRID_H
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include<vector> // for vector 
#include<algorithm> // for copy() and assign() 
#include<iterator> // for back_inserter 
#include <Misc.hpp>

using namespace std;

namespace grid
{

#define DEFAULT_PAGE_SIZE 4096


/**
 * @brief Class GridException Inheritance of class exception
 * @see RTMException.hpp
 */
class GridException : public exception
{
private:
    string what_msg; ///< Message for exception
public:
    /**
     * @brief Construct a new Grid Exception object
     * @param _what_msg 
     */
    GridException(string &_what_msg)
    {
        what_msg = _what_msg;
    }

    GridException(GridException &e)
    {
        what_msg = e.getWhatMsg();
    }
    ~GridException(){

    }
    GridException& operator=(GridException &e){
        what_msg = e.getWhatMsg();
        return *this;
    }
    string & getWhatMsg(){
        return what_msg;
    }

    /**
     * @brief Exception specifications
     * @return const char* string with explanatory information
     */    
    virtual const char *what() const throw()
    {
        return what_msg.c_str();
    }
 
};

/**
 * @brief Class Grid
 * @tparam GridData_t 
 */
template<typename GridData_t >
class Grid
{
protected:
    HostBuffer_t<GridData_t> *gridBuffer=nullptr;
    // GridData_t * grid = nullptr;
public:
    size_t gridSize  = 0; /// total size
    GridData_t MAXVAL;
    GridData_t MINVAL;
    /**
     * @brief Construct a new Grid object
     */
    Grid(){
        // grid = nullptr;
        gridSize = 0;
        gridBuffer = new HostBuffer_t<GridData_t>(1);
        MAXVAL = static_cast<GridData_t>(-10000000000.);
        MINVAL = static_cast<GridData_t>(10000000000.);
    }

    Grid(size_t _size)
    {
        gridBuffer = new HostBuffer_t<GridData_t>(1);
        Grid<GridData_t>::resizeGrid(_size);
        //grid = new GridData_t[_size];
        MAXVAL = static_cast<GridData_t>(-10000000000.);
        MINVAL = static_cast<GridData_t>(10000000000.);
    }

    /**
     * @brief Construct a new Grid object
     * @param g
     */
    Grid(const Grid<GridData_t> &g)
    {
        gridBuffer = new HostBuffer_t<GridData_t>(1);
        Grid<GridData_t>::resizeGrid(g.size());
        MAXVAL = g.MAXVAL;
        MINVAL = g.MINVAL;
        for (size_t i0 = 0; i0 < size; i0++)
        {
            (*this)[i0] = g[i0];
        }
    }
    /**
     * @brief Destroy the Grid object
     */
    ~Grid()
    {
        resizeGrid(0);
        delete gridBuffer;
    }
    /**
     * @brief  Assignment operator
     * @param g  reference to Grid<GridData_t>
     * @return Grid& 
     */
    Grid &operator=(const Grid<GridData_t> &g)
    {
        Grid<GridData_t>::resizeGrid(g.size());
        MAXVAL = g.MAXVAL;
        MINVAL = g.MINVAL;
        for (int i0 = 0; i0 < size; i0++)
        {
            (*this)[i0] = g[i0];
        }
        return *this;
    }
    GridData_t &operator[](int x) const
    {
        // return Grid<GridData_t>::grid[x];
        return Grid<GridData_t>::gridBuffer->at(x);
    }
    GridData_t *data() const
    {
        //return Grid<GridData_t>::grid;
        return Grid<GridData_t>::gridBuffer->data();
    }
    size_t size() const
    {
        return gridSize;
    }

    HostBuffer_t<GridData_t> & getGridBuffer(){
        return *Grid<GridData_t>::gridBuffer;
        // return gridBuffer;
    }
    /**
     * @brief Get the By Offset object
     * @param offset 
     * @return T& 
     */
    GridData_t &getByOffset(size_t offset) const
    {
        // return Grid<GridData_t>::grid[offset];
        return Grid<GridData_t>::gridBuffer->at(offset);
    }
    /**
     * @brief Set the By Offset object
     * @param offset 
     * @param val 
     */
    void setByOffset(size_t offset, GridData_t val)
    {
        // Grid<GridData_t>::grid[offset] = val;
        Grid<GridData_t>::gridBuffer->at(offset)  = val;
    }

    /**
     * @brief Function fill(): fills the entire grid data 
     *        with 'val'. Does not update device buffers.
     * @param val 
     */
    void fill(GridData_t val)
    {
        // std::fill(grid->begin(), grid->end(), val);
        size_t k;
        for (k=0; k<gridSize; k++){
            //Grid<GridData_t>::grid[k] = val;
            Grid<GridData_t>::gridBuffer->at(k) = val;
        }
    }
    /**
     * @brief Set the Max Min object
     */
    void setMaxMin()
    {
        size_t i0 = 0;
        for (i0 = 0; i0 < gridSize; i0++)
        {
            if (Grid<GridData_t>::gridBuffer->at(i0) >= MAXVAL)
            {
                MAXVAL = Grid<GridData_t>::gridBuffer->at(i0);
            }
            if (Grid<GridData_t>::gridBuffer->at(i0) <= MINVAL)
            {
                MINVAL = Grid<GridData_t>::gridBuffer->at(i0);
            }
        }
    }
    /**
     * @brief Function multiply()
     * @param m
     */
    void multiplyBy(GridData_t m)
    {
        size_t i0 = 0;
        for (i0 = 0; i0 < gridSize; i0++)
        {
            Grid<GridData_t>::gridBuffer->at(i0) *= m;
        }
    }

    void subtractBy(Grid<GridData_t> & _s){
        size_t i0 = 0;
        for (i0 = 0; i0 < gridSize; i0++)
        {
            Grid<GridData_t>::gridBuffer->at(i0) -= _s.getByOffset(i0);
        }
    }

    void power2(){
        size_t i0 = 0;
        for (i0 = 0; i0 < gridSize; i0++)
        {
            Grid<GridData_t>::gridBuffer->at(i0) *= Grid<GridData_t>::gridBuffer->at(i0);
        }
    }

    void resizeGrid(const size_t _size){
        Grid<GridData_t>::gridBuffer->resize(_size);
        // Grid<GridData_t>::grid = Grid<GridData_t>::gridBuffer->data();
        Grid<GridData_t>::gridSize = _size;
    }


    /**
     * @brief Function loadFromFile()
     * @param fname 
     */
    void loadFromFile(string &fname)
    {
        ifstream vfile(fname, ios::in | ios::binary);
        if (vfile.fail())
        {
            string msg = "Input grid file '" + fname + "' doesn't exist! Please check.";
            GridException ex(msg);
            throw ex; 
        }
        vfile.read(reinterpret_cast<char *>(Grid<GridData_t>::gridBuffer->data()), this->size() * sizeof(GridData_t));
        vfile.close();
    }
    /**
     * @brief Function loadFromFileAt()
     * @param fname 
     * @param cnt 
     * @param size 
     */
    void loadFromFileAt(string &fname, int cnt, size_t _size)
    {
        ifstream vfile(fname, ios::in | ios::binary);
        if (vfile.fail())
        {
            string msg = "Input grid file '" + fname + "' doesn't exist! Please check.";
            GridException ex(msg);
            throw ex; 
        }
        vfile.seekg(cnt * _size, ios::beg);
        vfile.read(reinterpret_cast<char *>(Grid<GridData_t>::gridBuffer->data()), this->size() * sizeof(GridData_t));
        vfile.close();
    }
    /**
     * @brief Function saveToFile()
     * @param fname 
     */
    void saveToFile(string &fname)
    {
        ofstream vfile(fname, ios::out | ios::binary);
        if (vfile.fail())
        {
            string msg = "Output grid file '" + fname + "' doesn't exist! Please check.";
            GridException ex(msg);
            throw ex; 
        }
        //vfile.write(reinterpret_cast<char *>(Grid<GridData_t>::gridBuffer->data()), this->size() * sizeof(GridData_t));
        vfile.write(reinterpret_cast<char *>(Grid<GridData_t>::gridBuffer->data()), this->size() * sizeof(GridData_t));
        vfile.close();
    }
    /**
     * @brief Function appendTofile()
     * @param fname  
     */
    void appendTofile(string &fname)
    {
        ofstream vfile(fname, ios::out | ios::binary | ios::app);
        if (vfile.fail())
        {
            string msg = "Output grid file '" + fname + "' doesn't exist! Please check.";
            GridException ex(msg);
            throw ex; 
        }
        vfile.write(reinterpret_cast<char *>(Grid<GridData_t>::gridBuffer->data()), gridSize * sizeof(GridData_t));
        vfile.close();
    }
    /**
     * @brief Function initGrid()
     * @param dims 
     * @param len 
     */
    void initGrid(size_t _size){
        resizeGrid(_size);
    }
    /**
     * @brief Function destroyGrid()
     */
    void destroyGrid(){

        Grid<GridData_t>::gridBuffer->resize(0);
        gridSize = 0;
    }

    /**
     * @brief Function destroyGrid()
     */
    void destroyGrid(GridData_t * l_grid){
        l_grid->resize(0);
    }

    void copyData(Grid<GridData_t> &_from){
        // Copying vector by copy function 
        Grid<GridData_t>::getGridBuffer().assign(_from.getGridBuffer().begin(), _from.getGridBuffer().end());
    }
};
} // namespace grid

#endif