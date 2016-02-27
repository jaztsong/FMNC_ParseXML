#include<stdio.h>
#include<sys/stat.h>
#include <stdlib.h>
#include<iostream>
#include<string.h>
#include<dirent.h>
#include<vector>
#include <omp.h>
#include "fmnc_parser.h"

using namespace std;
vector<string> listFile(const char* );
int getFileSize(const char* );
//class fmnc_parser;
int main(int argc, char* argv[]){
        vector<string> flist;
        if(argc > 2){
                cerr<<"Input param error: please input dir path."<<endl;
        }
        else if(argc == 1 ){
                string defaultPath("./");
                flist=listFile(defaultPath.c_str());
        }
        else{
                flist=listFile(argv[1]);
        }
        omp_set_num_threads(8);
#pragma omp parallel
        {
#pragma omp for
                for(uint32_t i=0;i<flist.size();i++){
                        fmnc_parser parser(flist[i]);
                        parser.dump_str();
                }
        }


        return 0;

}

int getFileSize(const char* fileName)
{
        struct stat statbuf;

        if (stat(fileName, &statbuf) == -1) {
                /* check the value of errno */
                cerr<<"*Error: Wrong file "<<fileName<<endl;
                return MAX_FILE_SIZE + 1;
        }
        return  (intmax_t) statbuf.st_size;

}
vector<string> listFile(const char* dir){
        string   str;
        DIR *pDIR;
        struct dirent *entry;
        vector<string> flist;
        cout<<"Opening dir: "<<dir<<endl;
        if(( pDIR=opendir(dir) ) != NULL){
                while(( entry = readdir(pDIR) )!= NULL){
                        string fn=entry->d_name;
                        size_t found = fn.find_last_of(".");
                        string ext =fn.substr(found+1);
                        if( strcmp(ext.c_str(), "xml") == 0 &&
                                        getFileSize((dir+fn).c_str()) < MAX_FILE_SIZE){
                                Debug(dir+fn<<": the File Size is "<<getFileSize((dir+fn).c_str()));
                                flist.push_back(dir+fn);
                        }

                }
                return flist; 
                closedir(pDIR);

        }
        else{
                cerr<<"***Open Dir error!!"<<endl;
        }
        return vector<string>();

}
