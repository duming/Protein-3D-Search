#include "Utility.hpp"


std::vector<std::string>
Utility:: split(std::string &s, char delim)
{
    std::vector<std::string> tokens;
    split(tokens, s, delim);
    return tokens;
}    



void 
Utility:: split(std::vector<std::string> &tokens, std::string &s, char delim )
{
    tokens.resize(0);
    std::stringstream ss(s); 
    std::string item;
    while (getline(ss, item, delim)) 
    {
        if(delim == ' ' && item.empty())
            continue;
        tokens.push_back(item);
    }
}




void
Utility:: split(std::vector<std::string> &tokens, const char* s, char delim)
{
    tokens.resize(0);
    do
    {
        if(' ' ==delim)
            while( ' ' == *s)
                s++;
        const char*begin = s;
        while(*s != delim && *s)
            s++;
        tokens.push_back(std::string(begin, s));
    }while(0 != *s++);

}



bool 
Utility:: listFiles(std::string dir, std::vector<std::string> &files, bool isFile)
{
	DIR *dp;
	struct dirent *dirp;
	if((dp  = opendir(dir.c_str())) == NULL)
    {
        std::cout << "Can't open " << dir << std::endl;
	    return false;
	}
	
	while ((dirp = readdir(dp)) != NULL) 
    {
        if((isFile && dirp->d_type == DT_REG)
            || (!isFile && dirp->d_type == DT_DIR)    
           )
	        files.push_back(std::string(dirp->d_name));
	}
	closedir(dp);
	return true;
}



////////////////////////////////////
//    readwrite buffer
///////////////////////////////////


template<class T>
const T& ReadWriteBuffer<T>:: startRead()
{
    full.P();
    pthread_mutex_lock(&rw_mutex);
    int ret = outptr;
    outptr = nextPos(outptr);
    return base[ret];
}




template<class T>
void ReadWriteBuffer<T>:: endRead()
{
    pthread_mutex_unlock(&rw_mutex);
    empty.V();
}


template<class T>
T& ReadWriteBuffer<T>:: startWrite()
{   
    empty.P();
    pthread_mutex_lock(&rw_mutex);
    int ret = inptr;
    inptr = nextPos(inptr);
    return base[ret];
    
}

template<class T>
void ReadWriteBuffer<T>:: endWrite()
{
    pthread_mutex_unlock(&rw_mutex);
    full.V();
}
