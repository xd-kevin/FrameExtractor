//
// Created by shawntang(唐兴) on 2019-02-25.
//

#ifndef FRAMEEXTRACTOR_FILE_HELPER_H
#define FRAMEEXTRACTOR_FILE_HELPER_H

#include <string>
#include <stdio.h>
#include <algorithm>


namespace extractor{

    struct FileParts
    {
        std::string path;
        std::string name;
        std::string ext;
    };

    static inline FileParts fileparts(std::string filename)
    {
        int idx0 = filename.rfind("/");
        int idx1 = filename.rfind(".");

        if( idx1 == (int) std::string::npos )
            idx1 = filename.length();

        FileParts fp;
        fp.path = filename.substr(0,idx0+1);
        fp.name = filename.substr(idx0+1,idx1-idx0-1);
        fp.ext  = filename.substr(idx1);

        return fp;
    };

    static inline std::string get_filename(std::string filepath){
        extractor::FileParts fp = extractor::fileparts(filepath);
        std::string filename = fp.name;
        replace(filename.begin(),filename.end(),' ','_');
        return filename;
    };
    static inline bool file_exsits(const std::string& name){
        if(FILE *file = fopen(name.c_str(),"r")){
            fclose(file);
            return true;
        }else{
            return false;
        }
    };
}
#endif //FRAMEEXTRACTOR_FILE_HELPER_H
