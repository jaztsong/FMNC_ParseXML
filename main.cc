#include<stdio.h>
#include<cstdlib>
#include<iostream>
#include<string.h>
#include<dirent.h>
#include<vector>
#include "fmnc_parser.h"

using namespace std;
vector<string> listFile(const char* );
class fmnc_parser;
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
        for(vector<string>::iterator it=flist.begin();it != flist.end();++it){
                fmnc_parser parser(*it);
                parser.dump_str();
        }

        return 0;

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
                        if( strcmp(ext.c_str(), "xml") == 0 ){
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
