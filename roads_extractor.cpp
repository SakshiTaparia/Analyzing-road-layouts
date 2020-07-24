#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <cstring>
#include <algorithm>
#include <stack>
#include <iomanip>
#include <fstream>
#include <queue>
#include <deque>
#include <stdio.h>

using namespace std;
#include<chrono>
using namespace std::chrono;
#include <random>
const long long max_nodes=4000000;
const long long max_ways=4000000;
const long long log_limit =25;
const long max_intersections=250;
const long double INF = 1e10;
const long long walk_ratio_index=50;
const long long epochs = 100;
string ignore_leading_space(string &s)
{
    string ret="";
    int i=0;
    while(i<s.size() and (s[i]==' ' or s[i]=='\t'))
    {
        i++;
    }
    while(i<s.size())
    {
        ret+=s[i];
        i++;
    }
    return ret;
}
vector<string> split_string(string s,char c)
{
    vector<string> list_of_str;
    int i=0,j=0;
    while(i<s.size())
    {
        string temp="";
        j=i;
        while(j<s.size() and (s[j])!=c and ((s[j])!='/') and  ((s[j])!='>'))
        {
            temp+=s[j];
            j++;
        }
        list_of_str.push_back(temp);
        i=j+1;
    }

    return list_of_str;
}
// The Knuth-Morris-Pratt algorithm for pattern matching
long long kmp(string s, string z)
{
    long long P=z.size();
    z+='#';
    z+=s;
    long long n=z.size();
    vector<long long> F(n);
    F[0]=0;
    for(int i=1;i<n;i++)
    {
        int j=F[i-1];
        while(j>0 and z[i]!=z[j])
            j=F[j-1];
        if(z[i]==z[j]) j++;
        F[i]=j;
    }
    long long ret=0;
    for(int i=(P+1);i<n;i++)
    {
        if(F[i]==P)
        {
            ret++;
        }
    }
    return ret;
}


string s;

//DATABASE
vector< map<string,string> > node_data;
vector< map<string,string> > way_data;

vector< vector< map<string,string> > > way_nodes,way_tags;
bool is_roundabout[max_ways],insignificant_node[max_nodes];
map<long long,long long> node_id_map,way_id_map;

//GRAPH


vector<long long> adjacency_list[max_nodes];
map<pair<long long,long long>,long double> edge_list;            //edge_list is a map which counts a particular edge's occurrences, we expect it to be 1 for the edges which
                                                     //are present in the graph and 0 otherwise
vector<pair<long long,long long>> list_of_edges;

long long degree[max_nodes];                         //stores the degrees of all the nodes


long long cnt_round;
void process_roundabouts()
{
    for(int i=0;i<way_tags.size();i++)
    {
        for(auto list_of_maps: way_tags[i])
        {
            for(auto key_value: list_of_maps)
            {
                if(kmp(key_value.first,"roundabout") or kmp(key_value.second,"roundabout"))
                {
                    is_roundabout[i]=true;
                    cnt_round++;
                    break;
                }
            }
            if(is_roundabout[i]) break;
        }
        if(is_roundabout[i])
        {
            string ID_string = way_nodes[i][0]["ref"];
            //cout<<ID_string<<"\n";
            ID_string=ID_string.substr(1,ID_string.size()-2);
            long long ID = stoll(ID_string);
            long long U = node_id_map[ID];
            for(int j=1;j<way_nodes[i].size()-1;j++)
            {
                string ID_string = way_nodes[i][j]["ref"];
                ID_string=ID_string.substr(1,ID_string.size()-2);
                long long ID = stoll(ID_string);
                long long V = node_id_map[ID];
                insignificant_node[V]=true;
                node_id_map[ID]=U;
            }
        }
    }
    cout<<cnt_round<<" roundabouts\n";
    //way_nodes.size()<<" "<<way_tags.size()<<"\n";

}
vector<string> ways_data,nodes_data;
vector<string> buffer;
int find_year(string &s){
    int y=0;
    for(auto u:s){
        y*=10;
        y+=(u-'0');
    }
    return y;
}
long long temp_count=0;
int revisions[5000];
void extract_roads()
{
	for(auto u:buffer)
	{
		cout<<u<<'\n';
	}
    for(int i=0;i<node_data.size();i++)
    {
        string year, timestamp = node_data[i]["timestamp"];
        year = timestamp.substr(1,4);
        revisions[find_year(year)]++;
        if(degree[i])
        {
            cout<<nodes_data[i]<<'\n';
        }
    }
    for(auto u:ways_data)
    {
        cout<<u<<'\n';
    }
	cout<<"</osm>";
}




void add_edges()
{
    int iter=0;
    for(auto way:way_nodes)
    {
        if(is_roundabout[iter++])
            continue;
        for(int i=0;i<way.size()-1;i++)
        {
            string ID_string1 = way[i]["ref"];
            string ID_string2 = way[i+1]["ref"];
            ID_string1=ID_string1.substr(1,ID_string1.size()-2);
            ID_string2=ID_string2.substr(1,ID_string2.size()-2);
            long long ID1 = stoll(ID_string1), ID2= stoll(ID_string2);
            long long U=node_id_map[ID1],V=node_id_map[ID2];
            pair<long long,long long> E1=make_pair(U,V),E2=make_pair(V,U);

            if((edge_list.find(E1)==edge_list.end()) and (edge_list.find(E2)==edge_list.end()) and (U!=V))
            {

                adjacency_list[U].push_back(V);
                adjacency_list[V].push_back(U);

                degree[U]++;
                degree[V]++;
                edge_list[E1]=1;
                edge_list[E2]=1;
                list_of_edges.push_back(E1);
                list_of_edges.push_back(E2);
            }
        }
    }

}

void way(ifstream & inp_handle)
{
    bool is_road = false;
    string temp_string = s;
    vector<string> way_strings = split_string(s,' ');
    map<string,string> way_data_map;
    way_data_map.clear();
    for(int i=1;i<way_strings.size()-1;i++)                        //first and last elements are not data elements
    {
        vector<string> temp = split_string(way_strings[i],'=');
        if(temp.size()>=2)
        way_data_map[temp[0]]=temp[1];

    }


    // store node data for every way
    vector<map<string,string>> cur_way_node_data,cur_way_tag_data;
    while(getline(inp_handle,s))
    {
        s=ignore_leading_space(s);
        temp_string+='\n';
        temp_string+=s;
        if(s.substr(1,2) =="nd")
        {
            vector<string> cur_way_node_strings = split_string(s,' ');
            map<string,string> cur_way_node_data_map;
            cur_way_node_data_map.clear();
            for(int i=1;i<cur_way_node_strings.size()-1;i++)             //first and last elements are not data elements
            {
                vector<string> temp = split_string(cur_way_node_strings[i],'=');
                if(temp.size()>=2)
                cur_way_node_data_map[temp[0]]=temp[1];
            }
            cur_way_node_data.push_back(cur_way_node_data_map);
        }
        else if(s.substr(1,3) =="tag")
        {

            vector<string> cur_way_tag_strings = split_string(s,' ');
            map<string,string> cur_way_tag_data_map;
            cur_way_tag_data_map.clear();
            for(int i=1;i<cur_way_tag_strings.size()-1;i++)             //first and last elements are not data elements
            {
                vector<string> temp = split_string(cur_way_tag_strings[i],'=');
                if(temp.size()>=2)
                {
                    cur_way_tag_data_map[temp[0]]=temp[1];

                    if(kmp(temp[1],"highway"))
                        is_road=true;
                }

            }

            cur_way_tag_data.push_back(cur_way_tag_data_map);
        }
        else
        {
            break;
        }

    }

    if(is_road)
    {
        way_nodes.push_back(cur_way_node_data);
        way_tags.push_back(cur_way_tag_data);
        ways_data.push_back(temp_string);
        string ID_string =way_data_map["id"];
        ID_string=ID_string.substr(1,ID_string.size()-2);
        long long ID = stoll(ID_string);

        way_id_map[ID]=(long long)way_data.size();                   //for almost O(1) access
        way_data.push_back(way_data_map);
    }
}
void node()
{
    if(s[s.size()-1]=='>' and s[s.size()-2]!='/')
    {
    	s+= "\n</node>";
    }
    nodes_data.push_back(s);

    vector<string> node_strings = split_string(s,' ');

    map<string,string> node_data_map;
    node_data_map.clear();
    //cout<<s<<endl;
    for(int i=1;i<node_strings.size()-1;i++)                      //first and last elements are ignored as they are node data elements
    {
        vector<string> temp = split_string(node_strings[i],'=');
        if(temp.size()>=2)
            node_data_map[temp[0]]=temp[1];
    }

    string ID_string =node_data_map["id"];
    ID_string=ID_string.substr(1,ID_string.size()-2);
    long long ID = stoll(ID_string);
    //cout<<ID<<endl;
    node_id_map[ID]=(long long)node_data.size();                 //for almost O(1) access
    node_data.push_back(node_data_map);
}

/**
Some useful conversions from string
stold -> string to long double
stoll -> string to long long
**/

int main()
{

    int no_of_files=1;
    //freopen("tess.txt","w",stdout);

        //auto start = high_resolution_clock::now();
        //cout<<fixed;
        //cout<<setprecision(4);


		node_data.clear();
        way_data.clear();
        way_nodes.clear();
        way_tags.clear();
        memset(is_roundabout,false,sizeof(is_roundabout));
        memset(insignificant_node,false,sizeof(insignificant_node));
        node_id_map.clear();
        way_id_map.clear();

        edge_list.clear();
        list_of_edges.clear();
        memset(degree,0LL,sizeof(degree));
        nodes_data.clear();
        ways_data.clear(); 
		buffer.clear();
        string inp_file ="",out_file ="";
            inp_file+="kolkata";
            // inp_file+=to_string(i);
            inp_file+=".osm";
            out_file+="kolkata_processed";
            // out_file+=to_string(i);
            out_file+=".osm";
            ifstream inp_handle;
            inp_handle.open(inp_file.c_str());
            //freopen(inp_file.c_str(),"r",stdin);
            freopen(out_file.c_str(),"w",stdout);
            //freopen("tess.txt","w",stdout);
		bool node_encountered=false;
        while(getline(inp_handle,s))
        {
                s=ignore_leading_space(s);
                //cout<<s<<endl;

                if(s.size()==0) continue;
                if(s.substr(1,4) =="node")
				{
					if(s[s.size()-1]=='>')
	                {
	                    s = s.substr(0,s.size()-1);
	                    if(s[s.size()-1]=='/')
	                    {
	                        s = s.substr(0,s.size()-1);
							s+=' ';
							s+='/';
							s+='>';
	                    }
						else
						{
			                s+=' ';
			                s+='>';
						}
	                }
                    node();
					node_encountered=true;
				}
                else if(s.substr(1,3) =="way")
                {
                	if(s[s.size()-1]=='>')
	                {
	                    s = s.substr(0,s.size()-1);
	                    if(s[s.size()-1]=='/')
	                    {
	                        s = s.substr(0,s.size()-1);
							s+=' ';
							s+='/';
							s+='>';
	                    }
						else
						{
			                s+=' ';
			                s+='>';
						}
	                }
                    way(inp_handle);
                }
                else
				{
                    if(!node_encountered)
					{
						buffer.push_back(s);
					}
				}
        }
        //process_roundabouts();
        for(int i=0;i<=node_data.size();i++)
        {
            adjacency_list[i].clear();
            //adj_dist[i].clear();
        }
        add_edges();
        extract_roads();

        freopen("mumbai_timestamps.csv","w",stdout);
        cout<<"Year,Revisions\n";
        for(int i=0;i<5000;i++){
            if(revisions[i]){
                cout<<i<<","<<revisions[i]<<'\n';
            }
        }
        cerr<<"done"<<endl;
    //identify_blocks();
 //    walkability_ratio();
 //    total_road_length();
    // find_blocks();
 //    auto end = high_resolution_clock::now();
 //    auto time_taken = duration_cast<microseconds>(end-start);
 //    cout<<"\nTime taken: "<<time_taken.count()<<" microseconds\n";

}
