#ifndef RTM_EXCEPTION_H
#define RTM_EXCEPTION_H

#include <cstdlib>
#include <cstddef>    // nullptr_t, ptrdiff_t, size_t
#include <string>     // string, stoi, to_string
#include <iostream>
#include <exception>


using namespace std; 
/**
 * @class Class RTMException Inheritance of class exception
 */
class RTMException : public exception
{
private:
    string what_msg; ///< Message for exception

public:
    /**
     * @brief Construct a new RTMException object
     * @param _what_msg 
     */
    RTMException(string &_what_msg)
    {
        what_msg = _what_msg;
    }
    RTMException(RTMException &e)
    {
        what_msg = e.getWhatMsg();
    }
    ~RTMException(){

    }
    RTMException& operator=(RTMException &e){
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


#endif