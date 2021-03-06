#include "fmnc_parser.h"

double average(vector<double> v)
{      double sum=0;
        for(uint32_t i=0;i<v.size();i++)
                sum+=v[i];
        return sum/v.size();
}
double fixTime(string s)
{
        vector<string> strs;
        boost::split(strs,s,boost::is_any_of("."));
        return atof(strs[0].c_str())+atof(strs[1].c_str())/1e6;
}
template <typename T>
std::string to_string(T value)
{
        std::ostringstream os;
        os << value;
        return os.str();
}
double pearson_correlation(vector<double> x,vector<double> y)
{
        double x_mean=average(x),y_mean=average(y);
        double upper_part=0,x_2_sum=0,y_2_sum=0;
        for(uint32_t i=0;i<x.size();i++){
                    upper_part +=( x[i]-x_mean )*(y[i] - y_mean);
                    x_2_sum += pow(x[i]- x_mean,2);
                    y_2_sum += pow(y[i]- y_mean,2);
        }
        return upper_part/sqrt(x_2_sum * y_2_sum);
}
/////////////////////////////////////////////////////////////////////////
/////////////////////UNAGGRE_CHUNK//////////////////////////////////
////////////////////////////////////////////////////////////////////////
void unaggre_chunk::setDataSet(fmnc_measurer_set* ms)
{
        if(strcmp(( ms->getLabel() ).c_str(),"Receive") == 0){
                mRcvdSet=ms;
        }
        else if( strcmp((ms->getLabel() ).c_str(),"Send") == 0){
                mSentSet=ms;
        }
}
void unaggre_chunk::set_length(uint16_t l)
{
        length = l;
}
uint16_t unaggre_chunk::get_length()
{
        return length;
}
uint32_t unaggre_chunk::get_start()
{
        return start_index;
}

void unaggre_chunk::set_j_end(uint32_t st)
{
        j_end = st;
}
void unaggre_chunk::set_j_start(uint32_t st)
{
        j_start = st;
}
void unaggre_chunk::find_jumbo()
{
        //Find Start
        uint32_t st=start_index,end = start_index+length-1;
        //default
        set_j_start(st);

        for(uint32_t offset=0;offset<(*mInterACK).size();offset++){
                Debug("looking start "<<st<<" "<<offset<<" "<<(*mInterACK).size());
                if(st < (*mInterACK).size() && st < mLab){
                        Debug(" + offset"<<st + offset);
                        if( st + offset <(*mInterACK).size() &&
                                         st + offset < mLab &&
                                        (*mInterACK)[st+offset] > AI_GAP_THRESHOLD ){
                                set_j_start(st + offset);
                                break;
                        }
                        Debug(" - offset "<<st - offset);
                        if(  st > offset &&
                                        (*mInterACK)[st-offset] > AI_GAP_THRESHOLD ){
                                set_j_start ( st - offset);
                                break;
                        }
                }
        }
        /* Debug("Find Start "<<st<<" "<<j_start); */
        //Find End
        /* for(int i = 0;i<(*mInterACK).size();i++) */
        /*         cout<<" InterACK "<<i<<" "<<(*mInterACK)[i]<<endl; */
        /* Debug("The End start from "<<end); */
        end = min((uint32_t)(*mInterACK).size(),end);
        end = min(mLab-2,end);
        //default
        set_j_end(end);
        for(uint32_t offset=0;offset+1<(*mInterACK).size() ;offset++){
                Debug("looking end "<<end<<" "<<offset<<" "<<(*mInterACK).size());
                if(end < (*mInterACK).size() - 1 && end < mLab ){
                        Debug(" + offset"<<end + offset);
                        if( end + offset <(*mInterACK).size() - 1 &&
                                        end + offset < mLab &&
                                        end + offset > j_start &&
                                        (*mInterACK)[end+offset] > AI_GAP_THRESHOLD ){
                                set_j_end(end + offset-1);
                                break;
                        }
                        Debug(" - offset "<<end - offset);
                        if( end  > j_start + offset &&
                                        (*mInterACK)[end-offset] > AI_GAP_THRESHOLD ){
                                set_j_end ( end - offset-1);
                                break;
                        }
                }
        }
        /* Debug("Find End "<<end<<" "<<j_end); */

}

void unaggre_chunk::calc_rates()
{
        if(j_start + 1 > j_end){
                cerr<<"ERROR:illegal unaggregated chunk with size "<<j_end -j_start<<endl;
                return ;
        }

        double sum_time=0;
        uint32_t sum_size=0;
        vector<fmnc_measurer_point*>* tmp = mSentSet->getData();
        /* cout<<"Start to calculating rates "<<j_start<<" "<<j_end<<endl; */
        Debug("Start to calculating rates "<<j_start<<" "<<j_end);

        
        //Send
        sum_time =1e6*( ( *tmp )[j_end+1]->get_time()-( *tmp )[j_start]->get_time() );
        
        for(uint32_t i = j_start;i<j_end+1;i++){
               sum_size += ( ( *tmp )[i]->get_size() )*8;
        }
        sent_rate = sum_size/sum_time;
        //receive
        tmp = mRcvdSet->getData();

        sum_time=0;
        sum_time =1e6*( ( *tmp )[j_end+1]->get_time()-( *tmp )[j_start]->get_time() );
        rcvd_rate = sum_size/sum_time;
        
        Debug(" rates: "<<start_index<<" rcvd: "<<rcvd_rate<<" send: "<<sent_rate);
        /* cout<<" rates: "<<start_index<<" rcvd: "<<rcvd_rate<<" send: "<<sent_rate<<endl; */
        
}
bool unaggre_chunk::decide_tag(float theta1,float theta2,uint8_t max)
{
        float theta = min(theta1*sent_rate,theta2);
        if(theta + rcvd_rate > sent_rate &&
                        sent_rate < max&&
                        rcvd_rate < sent_rate + 2){
                return true;
        }
        else
                return false;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////FMNC_MEASURER_POINT//////////////////////////////////
////////////////////////////////////////////////////////////////////////
//
fmnc_measurer_point::fmnc_measurer_point(double t,uint16_t s,uint32_t ts,uint32_t an)
{
        mTime=t;
        mSize=s;
        mTsVal=ts;
        mSN_AN=an;
        //TODO: Add more attribute later if needed For Lixing Thu 25 Feb 2016 01:43:03 PM EST.

}
double fmnc_measurer_point::get_time()
{
        return mTime;
}

uint16_t fmnc_measurer_point::get_size()
{
        return mSize;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////FMNC_MEASURER_SET//////////////////////////////////
////////////////////////////////////////////////////////////////////////
//
string fmnc_measurer_set::getLabel()
{
        return label;
}


void fmnc_measurer_set::add_item(fmnc_measurer_point* mp){
        mdata.push_back(mp);
}
vector<fmnc_measurer_point*>* fmnc_measurer_set::getData(){
        return &mdata;
}

///////////////////////////////////////////////////////////////////////////
/////////////////////FMNC_PARSER//////////////////////////////////
////////////////////////////////////////////////////////////////////////
fmnc_parser::fmnc_parser(string fn)
{
        mfilename=fn;
        mLab = DEFAULT_LAB;
        mRmax = DEFAULT_RMAX;
        mAB = 0.0;
        mRequestHelper.throughput = 0;
        mRequestHelper.accelerate = 0;
        init_times();
        init_AN_loger();
        


        start_parse();


}
void fmnc_parser::init_times()
{

        mTime_AB_end = 0;
        mTime_AB_start = 0;
        mTime_EI_end = 0;
        mTime_EI_start = 0;
        mTime_SYN_end = 0;
        mTime_SYN_start = 0;
}
void fmnc_parser::init_AN_loger()
{

        mAN_logger.nR1=1;
        mAN_logger.nR2=0;
        mAN_logger.nR3=0;
        mAN_logger.nR4=0;
        mAN_logger.nR5=0;
        mAN_logger.nE1=0;
        mAN_logger.nE2=0;
        mAN_logger.nOther=0;

}
void fmnc_parser::dump_str()
{
        
        string result("");
        //result += "filename=%s";
        //result += "time=%d";
        //result += "rtt=%.2f";
        //result += "pdr=%.2f";
        //result += "aggre=%.2f";
        //result += "ab=%.2f";
        //result += "cor=%.2f";
        //result += "ei=%.2f";
        //result += "ANs=%s";
        //result += "ab_time=%.2f";
        //result += "ei_time=%.2f";
        //result += "3way_time=%.2f";
        //result += "app=%s";
        //result += "id=%s";
        //result += "type=%s";
        //result += "ssid=%s";
        //result += "bssid=%s";
        //result += "rssi=%s";
        //result += "throughput=%u";
        //result += "latitude=%f.2f";
        //result += "longitude=%f.2f";
        //result += "accelerate=%f.2f\n";

        result += "%s,";
        result += "%d,";
        result += "%.2f,";
        result += "%.2f,";
        result += "%.2f,";
        result += "%.2f,";
        result += "%.2f,";
        result += "%.2f,";
        result += "%s,";
        result += "%.2f,";
        result += "%.2f,";
        result += "%.2f,";
        result += "%s,";
        result += "%s,";
        result += "%s,";
        result += "%s,";
        result += "%s,";
        result += "%s,";
        result += "%u,";
        result += "%s,";
        result += "%s,";
        result += "%.2f\n";

        printf(result.c_str(),get_filename().c_str(),getConnectionTime(),average(mRTT),
                        calc_packetloss(),calcAggregation(),getAB(),getCor(),getEI(),dump_AN_logger().c_str(),
                        getAB_duration(),getEI_duration(),get3way_duration(),
                        mRequestHelper.app.c_str(),mRequestHelper.id.c_str(),
                        mRequestHelper.type.c_str(),mRequestHelper.ssid.c_str(),
                        mRequestHelper.bssid.c_str(),mRequestHelper.rssi.c_str(),
                        mRequestHelper.throughput,mRequestHelper.latitude.c_str(),
                        mRequestHelper.longitude.c_str(),mRequestHelper.accelerate);
        dump_web();

}
struct tm * fmnc_parser::parse_time(uint64_t t)
{
        time_t tt=t;
        struct tm * ptm = localtime(&tt);
        return ptm;
}
void fmnc_parser::dump_web()
{
        struct tm *ptm=parse_time(getConnectionTime());
        string dir=WEB_RESULT_DIR;
        //Get time and date
        static  char date[20];
        sprintf(date,"%d-%.2d-%.2d",1900+ptm->tm_year,1+ptm->tm_mon,ptm->tm_mday);
        static  char time[20];
        sprintf(time,"%.2d:%.2d:%.2d",ptm->tm_hour,ptm->tm_min,ptm->tm_sec);

        dir += "/" +mRequestHelper.app+"-"+mRequestHelper.id;
        dir += "/" + string(date);
        Debug("The dump web dir "<<dir);
        boost::filesystem::path w_dir(dir);
        //Create the dir if not exist
        boost::filesystem::create_directories(w_dir);        
        string fn=dir+"/"+time+".html";
        ofstream myfile (fn.c_str());
        if (myfile.is_open())
        {
                myfile << prepare_web_content(ptm);
                myfile.close();
        }
        else 
                cerr << "**dump_web-Error:Unable to open file";

}
string fmnc_parser::prepare_web_content(struct tm * ptm)
{
        string ct = "";
        ct += "<!DOCTYPE html> \n\
                <html>\n\
                <head>\n\
                <title>\n\
                FMNC Test Result from";
        static  char time[20];
        static  char date[20];
        sprintf(date,"%d-%.2d-%.2d",1900+ptm->tm_year,1+ptm->tm_mon,ptm->tm_mday);
        sprintf(time,"%.2d:%.2d:%.2d",ptm->tm_hour,ptm->tm_min,ptm->tm_sec);
        ct += string(date)+" "+string(time);
        ct += " </title>";
	ct += "<script src=\"http://mbostock.github.com/d3/d3.v2.js\"></script>\n\
		<style>\n\
			path {\n\
				stroke: steelblue;\n\
				stroke-width: 1;\n\
				fill: none;\n\
			}\n\
			.axis {\n\
			  shape-rendering: crispEdges;\n\
			}\n\
			.x.axis line {\n\
			  stroke: lightgrey;\n\
			}\n\
			.x.axis .minor {\n\
			  stroke-opacity: .5;\n\
			}\n\
			.x.axis path {\n\
			  display: none;\n\
			}\n\
			.y.axis line, .y.axis path {\n\
			  fill: none;\n\
			  stroke: #000;\n\
			}\n\
		</style>";

        ct += " </head>\n\
                <body>\n\
                <h1>FMNC Test Result for "+string(date)+" "+string(time)+"</h1>\n\
                <table border=1>";

        string AB_result="";
        if(getAB()==12.0){
                AB_result = "> 10 Mb/s";
        }
        else if (getAB()==0.0){
                AB_result = "< 1 Mb/s";
        }
        else
                AB_result=to_string(getAB());

        ct +=  "<tr>\n\
                <td rowspan=\"6\"> User Profile</td>\n\
                <td> Client App </td>\n\
                <td>"+mRequestHelper.app+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> Client ID </td>\n\
                <td>"+mRequestHelper.id+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> SSID </td>\n\
                <td>"+mRequestHelper.ssid+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> BSSID </td>\n\
                <td>"+mRequestHelper.bssid+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> Longitude </td>\n\
                <td>"+mRequestHelper.longitude+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> Latitude </td>\n\
                <td>"+mRequestHelper.latitude+"</td>\n\
                </tr> ";

        ct +=  "<tr>\n\
                <td rowspan=\"6\"> FMNC Measurement</td>\n\
                <td> RTT (ms) </td>\n\
                <td>"+to_string(average(mRTT))+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> PDR </td>\n\
                <td>"+to_string(calc_packetloss())+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> Available Bandwidth (Mb/s) </td>\n\
                <td>"+AB_result+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> AB estimation Time (ms) </td>\n\
                <td>"+to_string(getAB_duration())+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> Elasticity Index </td>\n\
                <td>"+to_string(abs(getEI()-10))+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> Uplink Pearson Correlation </td>\n\
                <td>"+to_string(getCor())+"</td>\n\
                </tr>";
        ct +=  "<tr>\n\
                <td rowspan=\"3\">Achievable Throughput</td>\n\
                <td> Download Size (MB)</td>\n\
                <td>7.2</td>\n\
                </tr>\n\
                <tr>\n\
                <td> Download Time (ms)</td>\n\
                <td>"+to_string(mRequestHelper.throughput)+"</td>\n\
                </tr>\n\
                <tr>\n\
                <td> Average Throughput (Mb/s)</td>\n\
                <td>"+to_string(7.2*8*1000/mRequestHelper.throughput)+"</td>\n\
                </tr>";
        ct += "</table>\n\
		<h2>Analysis Graphs:</h2>";
	string data;
	stringstream ss;
	for(uint32_t i=0; i<mRTT.size();i++){
		//sprintf(buffer[i*6],"%.2f",mRTT[i]);
		ss<<mRTT[i]<<",";
	}
	data = ss.str();
	ct += prepare_D3JS("0","500",data,"RTT(ms)");

	ss.str("");
	ss.clear();
	for(uint32_t i=0; i<mInterACK.size();i++){
		//sprintf(buffer[i*6],"%.2f",mRTT[i]);
		ss<<mInterACK[i]<<",";
	}
	data = ss.str();
	ct += prepare_D3JS("1","850",data,"InterACK(ms)");

	ss.str("");
	ss.clear();
	for(uint32_t i=0; i<mlocalAggre.size();i++){
		//sprintf(buffer[i*6],"%.2f",mRTT[i]);
		ss<<mlocalAggre[i]<<",";
	}
	data = ss.str();
	ct += prepare_D3JS("1","850",data,"Local Aggregation Index");




        ct += "</body>\n\
               </html>";

        return ct;
                                    
}
string fmnc_parser::prepare_D3JS(string id, string pos, string data, string ylabel)
{
	string result="";
	result = "<div id=\"graph"+id+"\" class=\"aGraph\" style=\"position:absolute;top:"+pos+"px;left:0; float:left;\"></div>\n\
	<script>\n\
		var m = [80, 80, 80, 80]; \n\
		var w = 1000 - m[1] - m[3];\n\
		var h = 400 - m[0] - m[2]; \n\ 
		var data = [";
	result += data;

	result += "];\n\
		var x = d3.scale.linear().domain([0, data.length]).range([0, w]);\n\
		var y = d3.scale.linear().domain([0, d3.max(data)]).range([h, 0]);\n\
		var line = d3.svg.line()\n\
			.x(function(d,i) { \n\
				console.log('Plotting X value for data point: ' + d + ' using index: ' + i + ' to be at: ' + x(i) + ' using our xScale.');\n\
				return x(i); \n\
			})\n\
			.y(function(d) { \n\
				console.log('Plotting Y value for data point: ' + d + ' to be at: ' + y(d) + \" using our yScale.\");\n\
				return y(d); \n\
			})\n\
			var graph"+id+" = d3.select(\"#graph"+id+"\").append(\"svg:svg\")\n\
			      .attr(\"width\", w + m[1] + m[3])\n\
			      .attr(\"height\", h + m[0] + m[2])\n\
			    .append(\"svg:g\")\n\
			      .attr(\"transform\", \"translate(\" + m[3] + \",\" + m[0] + \")\");\n\
			var xAxis = d3.svg.axis().scale(x).tickSize(-h).tickSubdivide(true);\n\
			graph"+id+".append(\"svg:g\")\n\
			      .attr(\"class\", \"x axis\")\n\
			      .attr(\"transform\", \"translate(0,\" + h + \")\")\n\
			      .call(xAxis);\n\
			var yAxisLeft = d3.svg.axis().scale(y).ticks(4).orient(\"left\");\n\
			graph"+id+".append(\"svg:g\")\n\
			      .attr(\"class\", \"y axis\")\n\
			      .attr(\"transform\", \"translate(-25,0)\")\n\
			      .call(yAxisLeft);\n\
  			graph"+id+".append(\"svg:path\").attr(\"d\", line(data));\n\
			graph"+id+".append(\"text\")\n\
			.attr(\"text-anchor\", \"middle\")  // this makes it easy to centre the text as the transform is applied to the anchor\n\
			.attr(\"transform\", \"translate(\"+-60+ \",\"+(h/2)+\")rotate(-90)\")  // text is drawn off the screen top left, move down and out and rotate\n\
			.text(\""+ylabel+"\");\n\
			graph"+id+".append(\"text\")\n\
			.attr(\"text-anchor\", \"middle\")  // this makes it easy to centre the text as the transform is applied to the anchor\n\
			.attr(\"transform\", \"translate(\"+ (w/2) +\",\"+(h+(40))+\")\")  // centre below axis\n\
			.text(\"Packet No.\");\n\
	</script>";

        return result;
}

string fmnc_parser::get_filename()
{
        return mfilename;
}
double fmnc_parser::getAB_duration(){
        if(mTime_AB_end>0 && mTime_AB_start >0)
                //return millisecond
                return 1e3*(mTime_AB_end - mTime_AB_start);
        else
                return 1/0.0;
}

double fmnc_parser::getEI_duration(){
        if(mTime_EI_end>0 && mTime_EI_start >0)
                //return millisecond
                return 1e3*(mTime_EI_end - mTime_EI_start);
        else
                return 1/0.0;
}
double fmnc_parser::get3way_duration(){
        if(mTime_SYN_end>0 && mTime_SYN_start >0)
                //return millisecond
                return 1e3*(mTime_SYN_end - mTime_SYN_start);
        else
                return 1/0.0;
}
void fmnc_parser::start_parse()
{
        load_file(mfilename);
        if(parse_request()){
                calcRTT();
                calcInterACK();
                getLocalAggregation();
                find_unaggre_chunks();
                rate_analysis();
                detect_upbottleneck();
                calc_EI();
        }

}
string fmnc_parser::dump_AN_logger()
{
        char output[50];
        sprintf(output,"[%u:%u:%u:%u:%u:%u:%u:%u]",
                        mAN_logger.nR1,mAN_logger.nR2,
                        mAN_logger.nR3,mAN_logger.nR4,
                       mAN_logger.nR5,mAN_logger.nE1,
                      mAN_logger.nE2,mAN_logger.nOther);
        /* cout<<" "<<mAN_logger.nR1<<" "<<mAN_logger.nR2<<" "<<mAN_logger.nR3<<endl; */
        return string(output);
}
void fmnc_parser::update_AN_logger(string s)
{
        
                if(strcmp(s.c_str(),AN_R1) == 0)
                        mAN_logger.nR1++;
                else if(strcmp(s.c_str(),AN_R2) == 0)
                        mAN_logger.nR2++;
                else if(strcmp(s.c_str(),AN_R3) == 0)
                        mAN_logger.nR3++;
                else if(strcmp(s.c_str(),AN_R4) == 0)
                        mAN_logger.nR4++;
                else if(strcmp(s.c_str(),AN_R5) == 0)
                        mAN_logger.nR5++;
                else if(strcmp(s.c_str(),AN_E1) == 0)
                        mAN_logger.nE1++;
                else if(strcmp(s.c_str(),AN_E2) == 0)
                        mAN_logger.nE2++;
                else
                        mAN_logger.nOther++;
        
}
void fmnc_parser::load_file(string fn)
{
        pugi::xml_document doc;

        pugi::xml_parse_result result = doc.load_file(fn.c_str());

        if(result >0){
                Debug("Start to load file "<<fn.c_str());
                string time= doc.child("ConnectionTCPSlice").attribute("CreationTime").value();
                /* uint64_t tt = std::atoi(time.c_str()); */
                setRequest(doc.child("ConnectionTCPSlice").attribute("Request").value());
                setConnectionTime(std::atoi(time.c_str()));
                pugi::xml_node st = doc.child("ConnectionTCPSlice").child("MeasureRcvd");
                fmnc_measurer_set* ms;
                ms = new fmnc_measurer_set("Receive");
                fmnc_measurer_point* mp;
                uint16_t p_counter=0;
                for(pugi::xml_node_iterator it =st.begin();it != st.end();++it){
                        p_counter++;
                        if((strcmp( it->name(),"PktTCP") == 0 ) && 
                                        ( strcmp(it->attribute("Meta").value(),"") != 0  )) {
                                mp =new fmnc_measurer_point(fixTime(it->attribute("Time").value()),
                                                std::atoi(it->attribute("IPLength").value()),
                                                std::atoi(it->attribute("TsVal").value()),
                                                std::atoi(it->attribute("AN").value()));
                                ms->add_item(mp);

                                //Update the logger for the purpose of parsing loss later
                                update_AN_logger(it->attribute("AN").value());
                                /* if(( strcmp(it->attribute("AN").value(),LAST_AB_ACK_AN) == 0 ) && */
                                /*                 fixTime(it->attribute("Time").value())>mTime_AB_end) { */
                                /*         setAB_end(fixTime(it->attribute("Time").value())); */
                                /* } */
                                if(p_counter == 100) {
                                        setAB_end(fixTime(it->attribute("Time").value()));
                                }
                                else if( fixTime(it->attribute("Time").value())>mTime_EI_end) {
                                        setEI_end(fixTime(it->attribute("Time").value()));
                                }
                        }
                }
                setDataSet(ms);
                st = doc.child("ConnectionTCPSlice").child("MeasureSent");
                ms = new fmnc_measurer_set("Send");
                for(pugi::xml_node_iterator it =st.begin();it != st.end();++it){
                        if(strcmp( it->name(),"PktTCP") == 0  ) {
                                mp =new fmnc_measurer_point(fixTime(it->attribute("Time").value()),
                                                std::atoi(it->attribute("IPLength").value()),
                                                std::atoi(it->attribute("TsVal").value()),
                                                std::atoi(it->attribute("AN").value()));
                                ms->add_item(mp);
                                if( strcmp(it->attribute("SN").value(),FIRST_AB_SN) == 0 ) {
                                        setAB_start(fixTime(it->attribute("Time").value()));
                                }
                                else if(( strcmp(it->attribute("SN").value(),FIRST_EI_SN) == 0 )) {
                                        setEI_start(fixTime(it->attribute("Time").value()));
                                }

                        }
                }
                setDataSet(ms);
                //Capture the 3-way shake info
                st = doc.child("ConnectionTCPSlice").child("MeasureSetup");
                for(pugi::xml_node_iterator it =st.begin();it != st.end();++it){
                        if(strcmp( it->name(),"PktSYN") == 0  ) {
                                setSYN_start(fixTime(it->attribute("Time").value()));
                        }
                        else if((strcmp( it->name(),"PktTCP") == 0 ) && 
                                        ( strcmp(it->attribute("Meta").value(),"DataIn") == 0  )) {
                                setSYN_end(fixTime(it->attribute("Time").value()));
                        }
                }
                /* std::cout << "Load result: " << result << ",  Connection Request " << doc.child("ConnectionTCPSlice").attribute("Request").value() << std::endl; */
        }
}


fmnc_measurer_set* fmnc_parser::getDataSet(string s)
{
        if(strcmp(s.c_str(),"Receive") == 0)
                return mRcvdSet;
        else if(strcmp(s.c_str(),"Send") == 0)
                return mSentSet;
        else{
                cerr<<"*Error: getDataSet from parser: Wrong String"<<endl;
                return  0;
        }
}
void fmnc_parser::setDataSet(fmnc_measurer_set* ms)
{
        if(strcmp(( ms->getLabel() ).c_str(),"Receive") == 0){
                mRcvdSet=ms;
        }
        else if( strcmp((ms->getLabel() ).c_str(),"Send") == 0){
                mSentSet=ms;
        }
}
bool fmnc_parser::parse_request()
{
        Debug("Start to parse the request "<<mRequest);
        if(mRequest.find("/tests/train?0?SliceSize=") == std::string::npos){
                cerr<<"*ERROR: wrong request: "<<mRequest<<endl;
                return false;
        }
        vector<string> strs;
        boost::split(strs,mRequest,boost::is_any_of("?"));
        for(vector<string>::iterator it=strs.begin();it != strs.end();++it){
                /* Debug("Parseing request "<<*it); */
                if((*it).find("Length=") != std::string::npos)
                        mLab = atoi((*it).substr(7).c_str());
                else if((*it).find("Rmax=") != std::string::npos)
                        mRmax = atoi((*it).substr(5).c_str());
                else if((*it).find("app=") != std::string::npos)
                        mRequestHelper.app = (*it).substr(4);
                else if((*it).find("uuid=") != std::string::npos ||
                                (*it).find("imei=") != std::string::npos)
                        mRequestHelper.id = (*it).substr(5);
                else if((*it).find("type=") != std::string::npos)
                        mRequestHelper.type = (*it).substr(5);
                else if((*it).find("SSID=") != std::string::npos &&
                           (*it).find("SSID=") == (size_t) 0 )
                        mRequestHelper.ssid = (*it).substr(5);
                else if((*it).find("BSSID=") != std::string::npos &&
                           (*it).find("BSSID=") == (size_t) 0 )
                        mRequestHelper.bssid = (*it).substr(6);
                else if((*it).find("RSSI=") != std::string::npos)
                        mRequestHelper.rssi = (*it).substr(5);
                else if((*it).find("Throughput=") != std::string::npos)
                        mRequestHelper.throughput= atoi((*it).substr(11).c_str());
                else if((*it).find("Lat=") != std::string::npos)
                        mRequestHelper.latitude = (*it).substr(4).c_str();
                else if((*it).find("Long=") != std::string::npos)
                        mRequestHelper.longitude = (*it).substr(5).c_str();
                else if((*it).find("Acc=") != std::string::npos)
                        mRequestHelper.accelerate = atof((*it).substr(4).c_str());
                else if((*it).find("TPError=true") != std::string::npos)
                        mRequestHelper.throughput= 0;
                //TODO: The crude way to filter out throughput bad data need to be changed For Lixing Sat 05 Mar 2016 06:12:11 PM EST.
                
        }
        return true;

}
void fmnc_parser::setConnectionTime(uint64_t t)
{
        mConnectionTime = t;
}
uint64_t fmnc_parser::getConnectionTime()
{
        return mConnectionTime;
}
void fmnc_parser::calcRTT()
{
        Debug("Start to calculate RTT");
        vector<fmnc_measurer_point*>* sent = mSentSet->getData();
        vector<fmnc_measurer_point*>* rcvd = mRcvdSet->getData();
        for(uint32_t i=0;i+1<( *sent ).size() && i+1<( *rcvd ).size() ;i++){
                /* cout<<"Sent Time "<<( *sent )[i]->get_time()<<" Rcvd Time "<<( *rcvd )[i]->get_time( )<<endl; */
                mRTT.push_back(1e3*( ( *rcvd )[i]->get_time()- ( *sent )[i]->get_time( )));//Change to milisecond
        }
        Debug("Calculated the RTT with size "<<mRTT.size());
}
void fmnc_parser::calcInterACK()
{
        vector<fmnc_measurer_point*>* rcvd = mRcvdSet->getData();
        for(uint32_t i=0;i+2<( *rcvd ).size() ;++i){
                mInterACK.push_back(1e3*( ( *rcvd )[i+1]->get_time()- ( *rcvd )[i]->get_time( )));//Change to milisecond
        }
        Debug("Calculated the InterACK with size "<<mInterACK.size());
}


float fmnc_parser::calcAggregation()
{
        if(mInterACK.size()>0){
                int ag=0,sum=0;
                for(uint32_t i=0;( i+1<mInterACK.size() ) && ( i+1 < mLab );i++){
                        if(mInterACK[i]<float( AI_GAP_THRESHOLD ))
                                ag++;
                        sum++;
                }
                return float(ag)/float(sum);

        }
        else{
                cerr<<mfilename<<"calcAggregation:*Error: InterAck is empty"<<endl;
                return 1/0.0;
        }
}
void fmnc_parser::getLocalAggregation()
{
      
        if(mInterACK.size()>0){
                for(uint32_t i=0;mInterACK.size()>i + DEFAULT_LOCAL_WINDOW+1
                                &&i +DEFAULT_LOCAL_WINDOW+1< mLab;i++){
                        int ag=0;
                        for(uint32_t j=0;j<DEFAULT_LOCAL_WINDOW;j++){
                                if(mInterACK[i+j]<float( AI_GAP_THRESHOLD ))
                                        ag++;
                        }
                        mlocalAggre.push_back(float(ag)/float(DEFAULT_LOCAL_WINDOW));
                        /* cout<<i<<" "<<float(ag)/float(DEFAULT_LOCAL_WINDOW)<<endl; */
                }

        }
        else{
                cerr<<mfilename<<"getLocalAggregation:*Error: InterAck is empty"<<endl;
        }
        Debug("Calculated local Aggregation list with size "<<mlocalAggre.size());
    
}
void fmnc_parser::find_unaggre_chunks()
{
        uint32_t start=0,pre_unaggre=0,length=0;
        unaggre_chunk* uc;
        
        for(uint32_t i=0;i<mlocalAggre.size();i++)
        {
                if(mlocalAggre[i]<AI_THRESHOLD ){
                        //Find unaggregated part
                        /* cout<<"Unaggre: "<<i<<" "<<mlocalAggre[i]<<" "<<i-pre_unaggre<<endl; */
                        
                        if(i>AI_HOLE_LENGTH + pre_unaggre){
                                //If the unaggregation discontinued, then push what
                                //we have into the chunks
                                /* cout<<"Break"<<start<<" "<<length<<endl; */
                                
                                while(length > CHUNK_SIZE){
                                        uc = new unaggre_chunk(start);
                                        uc->set_length(CHUNK_SIZE);
                                        start += CHUNK_SIZE;
                                        length -= CHUNK_SIZE;
                                        m_unaggre_chunks.push_back(uc);
                                }
                                if(length > MIN_CHUNK_SIZE){
                                        uc = new unaggre_chunk(start);
                                        uc->set_length(length);
                                        m_unaggre_chunks.push_back(uc);
                                }
                                /* uc->setDataSet(getDataSet("Receive")); */
                                /* uc->setDataSet(getDataSet("Send")); */
                                length=1;
                                start=i;
                        } else{
                                //if the unaggregation part continues
                                if(length == 0)
                                        start = i;
                                length+=max(i-pre_unaggre,(uint32_t) 1);
                        }
                        pre_unaggre = i;
                }
                if(i+1 == mlocalAggre.size() && length > MIN_CHUNK_SIZE){
                        /* cout<<"Final Break"<<start<<" "<<length<<endl; */
                        while(length > CHUNK_SIZE){
                                uc = new unaggre_chunk(start);
                                uc->set_length(CHUNK_SIZE);
                                start += CHUNK_SIZE;
                                length -= CHUNK_SIZE;
                                m_unaggre_chunks.push_back(uc);
                        }
                        if(length > 0){
                                uc = new unaggre_chunk(start);
                                uc->set_length(length);
                                m_unaggre_chunks.push_back(uc);
                        }

                }

        }
}
void fmnc_parser::rate_analysis()
{
        Debug("Anaysis the rate:");
        uint32_t accu_false = 0;
        uint32_t end = 0;
        for(vector<unaggre_chunk*>::iterator it=m_unaggre_chunks.begin();it !=m_unaggre_chunks.end();++it){
                Debug(" Chunk "<<(*it)->get_start()<<" /w length=" <<(*it)->get_length()); 
                (*it)->setLab(get_mLab());
                (*it)->setDataSet(getDataSet("Receive"));
                (*it)->setDataSet(getDataSet("Send"));
                (*it)->setInterACK(&mInterACK);
                (*it)->find_jumbo();
                (*it)->calc_rates();
                if(!(*it)->decide_tag(THETA1,THETA2,mRmax)){
                        accu_false +=(*it)->get_length();
                        if(accu_false > ACCU_FALSE_LMT){
                                break;
                        }
                        Debug(" Bad "<<accu_false);
                }else{
                        mAB = (*it)->get_rcvd_rate();
                        end = (*it)->get_start()+(*it)->get_length()-1;
                        /* cout<<" Good "<<endl; */
                        Debug(" Good ");
                }
        }
        if(mAB > MAX_RATE_TARGET){
                mAB = 12;
                return;
        }
        //Check the end part
        unaggre_chunk uc = unaggre_chunk(end);
        uc.set_length(min(get_mLab(),(uint32_t)mInterACK.size()));//It does not matter since we have filter later
        uc.setLab(get_mLab());
        uc.setDataSet(getDataSet("Receive"));
        uc.setDataSet(getDataSet("Send"));
        uc.setInterACK(&mInterACK);
        Debug("Final Evaluate Chunk "<<uc.get_start()<<" /w length=" <<uc.get_length());
        uc.find_jumbo();
        uc.set_j_end((uint32_t)(min((uint32_t)mLab-2,(uint32_t)mInterACK.size())));
        uc.calc_rates();
        if(uc.decide_tag(THETA1,F_THETA2,mRmax)){
                Debug(" Good ");
                mAB = 12;
                return;
        }
        /* } */
        
}
void fmnc_parser::detect_upbottleneck()
{
        Debug("Detect Up:");
        std::vector<fmnc_measurer_point*>* rcvd_list = mRcvdSet->getData() ;
        //Find the common non zero segment 
        uint32_t t_len = min(mLab,(uint32_t)(*rcvd_list).size());
        std::vector<double> v_tsp,v_rcvd;
        /* cout<<"Start to Populate "<<t_len<<endl; */
        if(t_len < 1){
                cerr<<mfilename<<"*Error: Detect up find length < 1"<<endl;
                return;
        }
        
        Debug("The base for SENTACK is "<<float((*rcvd_list)[t_len-1]->getTsVal() - (*rcvd_list)[0]->getTsVal()));
        Debug("The base for RCVD is "<<float((*rcvd_list)[t_len-1]->get_time() - (*rcvd_list)[0]->get_time()));
        for(uint32_t i=0;i+2<t_len ;i++){
                v_tsp.push_back(((int32_t)( (*rcvd_list)[i+1]->getTsVal() ) - (int32_t)( (*rcvd_list)[i]->getTsVal() ))/
                                float((int32_t)( (*rcvd_list)[t_len-1]->getTsVal() ) -(int32_t) ( (*rcvd_list)[0]->getTsVal( ))));
                /* cout<<" Push timestamp "<<i<<" "<<((int32_t)( (*rcvd_list)[i+1]->getTsVal() ) - (int32_t)( (*rcvd_list)[i]->getTsVal() )) */
                /*         << " "<< (int32_t)( (*rcvd_list)[i]->getTsVal() ) <<endl; */
                v_rcvd.push_back(( (*rcvd_list)[i+1]->get_time()  - (*rcvd_list)[i]->get_time() )/
                                float( (*rcvd_list)[t_len-1]->get_time()  -  (*rcvd_list)[0]->get_time() ));
                /* cout<<" Push received "<<i<<" "<<( (*rcvd_list)[i+1]->get_time()  -  (*rcvd_list)[i]->get_time() ) */
                /*         << " "<<  (*rcvd_list)[i]->get_time()  <<endl; */
        }
        /* cout<<"Finish populate the lists "<<v_tsp.size()<<" "<<v_rcvd.size()<<endl; */
        
        mCor = pearson_correlation(v_tsp, v_rcvd);
        Debug("Correlation "<<mCor);
        
}
double fmnc_parser::calc_packetloss()
{
        std::vector<fmnc_measurer_point*>* rcvd_list = mRcvdSet->getData() ;
        return (uint32_t)(*rcvd_list).size()/float(mLab + EI_LENGTH);
        
}
void fmnc_parser::calc_EI()
{
        
        vector<double> vEI;
        if(mInterACK.size()<mLab){
                mEI = 1/0.0;
                return;
        }
        else{
                for(uint32_t i=mLab;i<mInterACK.size();i++){
                        if(mInterACK[i] > AI_GAP_THRESHOLD)
                                vEI.push_back(mInterACK[i]);
                        /* cout<<"EI: "<<i<<" "<<mInterACK[i]<<endl; */
                }
                mEI = average(vEI);
        }

}

