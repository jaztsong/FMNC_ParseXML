#include<stdio.h>
#include<cstdlib>
#include <stdint.h>
#include<iostream>
#include<typeinfo>
#include <string>
#include <string.h>
#include<sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include<vector>
#include <cmath>
#include "pugixml-1.7/src/pugixml.hpp"
using namespace std;

//Debug
//#define DEBUG

#ifdef DEBUG
#define Debug( x  ) std::cout << x <<endl;
#else
#define Debug( x  ) 
#endif

#define MAX_FILE_SIZE 80000
#define EI_LENGTH  40
#define DEFAULT_RMAX  15
#define MIN_CHUNK_SIZE  3
#define AI_HOLE_LENGTH  5
#define AI_GAP_THRESHOLD  0.4f
#define DEFAULT_LAB  100
#define DEFAULT_LOCAL_WINDOW 5
#define AI_THRESHOLD 0.8f
#define CHUNK_SIZE 10
#define ACCU_FALSE_LMT 19
#define THETA1 0.1f
#define THETA2 0.5
#define MAX_RATE_TARGET 10


class fmnc_measurer_point
{
        public:
                fmnc_measurer_point (double t,uint16_t s, uint32_t ts);
                double get_time();
                uint16_t get_size();
                uint32_t getTsVal(){return mTsVal;};

        private:
                double mTime;
                uint16_t mSize;
                uint32_t mTsVal;
                /* data */

};
class fmnc_measurer_set
{
        public:
                fmnc_measurer_set (string s):label(s){};
                void add_item(fmnc_measurer_point* mp);
                string getLabel();
                vector<fmnc_measurer_point*>* getData();


        private:
                string label;
                /* data */
                std::vector<fmnc_measurer_point*> mdata;

};
class unaggre_chunk
{
        public:
                unaggre_chunk (uint32_t st):start_index(st){};
                void set_length(uint16_t l);
                void calc_rates();
                void find_jumbo();
                uint32_t get_start();
                uint16_t get_length();
                void setDataSet(fmnc_measurer_set*);
                void setInterACK(vector<double>* m){mInterACK = m;};
                void setLab(uint32_t t){mLab = t;};
                bool decide_tag(float theta1,float theta2,uint8_t max);
                void set_j_end(uint32_t st);
                void set_j_start(uint32_t st);
                float get_rcvd_rate(){return rcvd_rate;};
        private:

                uint32_t mLab;
                uint32_t j_start;
                uint32_t j_end;
                uint32_t start_index;
                uint16_t length;
                float sent_rate;
                float rcvd_rate;
                fmnc_measurer_set* mSentSet;
                fmnc_measurer_set* mRcvdSet;
                std::vector<double>* mInterACK;

                /* data */

};

class fmnc_parser
{
        public:
                fmnc_parser(string fn);
                void dump_str();
                fmnc_measurer_set* getDataSet(string s);
                uint32_t get_mLab(){return mLab;};
        private:
                struct request_helper{
                        string app;
                        string id;
                        string type;
                        string ssid;
                        string bssid;
                        string rssi;
                        uint32_t throughput;
                        float latitude;
                        float longitude;
                        float accelerate;
                };
                request_helper mRequestHelper;
                string mRequest;
                uint32_t mLab;
                uint8_t mRmax;
                float mAB;
                double mCor;
                double mEI;
                std::vector<float> mlocalAggre;
                std::vector<unaggre_chunk*> m_unaggre_chunks;
                string mfilename;
                string get_filename();
                fmnc_measurer_set* mRcvdSet;
                fmnc_measurer_set* mSentSet;
                std::vector<double> mRTT;
                std::vector<double> mInterACK;



                void parse_request();
                void setRequest(string s){mRequest = s;};
                void calcRTT();
                double calc_packetloss();
                float calcAggregation(uint16_t st,uint16_t ed);
                void calcInterACK();
                float calcAggregation();
                void getLocalAggregation();
                void find_unaggre_chunks();
                void setDataSet(fmnc_measurer_set* ms);
                void rate_analysis();
                float getAB(){return mAB;};
                float getCor(){return mCor;};
                float getEI(){return mEI;};
                void detect_upbottleneck();
                void calc_EI();
                


                
                uint64_t mConnectionTime;
                void start_parse();
                void load_file(string fn);
                uint64_t getConnectionTime();
                void setConnectionTime(uint64_t t);



};
