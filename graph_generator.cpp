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
#include<set>
using namespace std;
#include<chrono>
using namespace std::chrono;
#include <random>
const long long max_nodes=800000;
const long long max_ways=100000;
const long long log_limit =25;
const long max_intersections=250;
const long double INF = 1e10;
const long long walk_ratio_index = 20;
const long long epochs = 10;
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
//stores the nodes(or ways) data for each of the grids
vector< map<string,string> > node_data[250000],node_store;
// vector< map<string,string> > way_data[250000];
map<pair<long long, long long>,long double> edge_list_all[250000];

long long rows,cols; // number of rows and columns in the rectangle   
// cols = (maxlat-minlat)*100.0

long double maxlat,minlat,maxlon,minlon;

bool usable_result[max_nodes];

vector< vector< map<string,string> > > way_nodes,way_tags;

// bool is_roundabout[max_ways],insignificant_node[max_nodes];
// roundabouts are being excluded from the current discussion now!


map<long long,long long> node_id_map,way_id_map;
map<long long, long long> current_map;      //maintains the map for the current grid


//GRAPH DATA

/**
This data will be used for every graph which gets created for every grid.
**/

vector<long long> adjacency_list[max_nodes];
vector<pair<long long,long double>> adj_dist[max_nodes];    //stores the adjacent vertices along with their corresponding edge weights
// map<pair<long long,long long>,long double> edge_list;            //edge_list is a map which counts a particular edge's occurrences, we expect it to be 1 for the edges which
                                                     //are present in the graph and 0 otherwise
vector<pair<long long,long long>> list_of_edges;

long long degree[max_nodes];                         //stores the degrees of all the nodes

/** A function to find the distance between two points given their latitude and longitude. We use the Haverline Formula **/


//RESULTS

vector<long long> res_three_way,res_four_way,res_more_than_four_way;
vector<long double> res_road_length, res_walkability_ratio;








long double find_distance(long double lat1,long double lon1, long double lat2,long double lon2)
{
    //all calculations are in SI units
    long double del_lat,del_lon,A,C,D;
    long double pi = acos(-1.0), rad_earth = 6378100.0;
    lat1*=(pi/180.0);
    lat2*=(pi/180.0);
    lon1*=(pi/180.0);
    lon2*=(pi/180.0);
    del_lat = lat2-lat1;
    del_lon = lon2-lon1;
    A = 0.5*(1-cos(del_lat)) + 0.5*cos(lat1)*cos(lat2)*(1-cos(del_lon));
    C = 2*atan2(sqrtl(A),sqrtl(1-A));
    D = rad_earth*C;
    return D;
}

void floyd_warshall(long double dist[][walk_ratio_index+1])
{
    long double eps = 0.0000000001;
    long long W=walk_ratio_index;
    for(int k=0;k<W;k++)
    {
        for(int i=0;i<W;i++)
        {
            for(int j=0;j<W;j++)
            {
                if(dist[i][k]+dist[k][j]<dist[i][j]-eps)
                {
                    dist[i][j]=dist[i][k]+dist[k][j];
                }
            }
        }
    }
}

void dijkstra(long long N, long long s, vector<long double> & d, vector<long long> & p, long long pos)
{
    // cerr<<"enter dijkstra"<<endl;
    // long long n = node_data[N].size();
    d.assign(N, INF);
    p.assign(N, -1);

    d[s] = 0;
    priority_queue<pair<long double,long long>> q;
    q.push({0, s});
    while (!q.empty())
    {
        long long v = q.top().second;
        long double d_v = -(q.top().first);
        q.pop();
        if (abs(d_v-d[v]) > 0.00000001)
            continue;
        if(v==pos) break;       // if we have found the shortest distance for the needed vertex, we can break;
        for (auto edge : adj_dist[v])
        {
            long long to = edge.first;
            long double len = edge.second;

            if (d[v] + len < d[to])
            {
                d[to] = d[v] + len;
                p[to] = v;
                q.push({-d[to], to});
            }
        }
    }
    // cerr<<"exit dijkstra"<<endl;

}

// finds the lat value given the vertex number
long double find_lat(long long v)
{
    string lat_val = node_store[v]["lat"];
    long double lat = stold(lat_val.substr(1,lat_val.size()-2));
    return lat;
}

// finds the lon value given the vertex number
long double find_lon(long long v)
{
    string lon_val = node_store[v]["lon"];
    long double lon = stold(lon_val.substr(1,lon_val.size()-2));
    return lon;
}

/** A straightforward depth first search **/

void dfs(vector<long long> adj[],bool vis[],long long u, long long p,map<pair<long long,long long>,long long>& edge_list,vector<long long> tree[])
{
    //cout<<node_data[u]["id"]<<endl;
    vis[u]=true;
    for(long long x:adj[u])
    {
        if(!vis[x])
        {
            tree[u].push_back(x);
            tree[x].push_back(u);
            vis[x]=true;
            dfs(adj,vis,x,u,edge_list,tree);
        }
        else
        {
            if(x!=p)
            {
                if(u>x)
                    edge_list[make_pair(u,x)]++;
                else
                    edge_list[make_pair(x,u)]++;
            }
        }
    }
}

long long time_in[max_nodes],time_out[max_nodes],up[max_nodes][log_limit+2],timer;

void lca_dfs(long long v,long long p,vector<long long> tree[],vector<long long>& list_of_vertices)
{
    list_of_vertices.push_back(v);
    time_in[v]=++timer;
    up[v][0]=p;
    for(int i=1;i<=log_limit;i++)
    {
        up[v][i] = up[up[v][i-1]][i-1];
    }
    for(long long u:tree[v])
    {
        if(u!=p)
        {
            lca_dfs(u,v,tree,list_of_vertices);
        }
    }
    time_out[v]=++timer;
}

bool is_ancestor(long long u,long long v)
{
    return ((time_in[u]<=time_in[v]) and (time_out[u]>=time_out[v]));
}

long long lca(long long u,long long v)
{
    if(is_ancestor(u,v)) return u;
    if(is_ancestor(v,u)) return v;
    for(int i=log_limit;i>=0;i--)
    {
        if(!is_ancestor(up[u][i],v))
            u = up[u][i];
    }
    return up[u][0];
}

/** The function to find the cycle splitter **/
void find_cycle_split(long long u,vector<long long> tree[],map<long long,long long> & node_in_cycle,vector<long long> & split,bool & done,map<long long,long long>& vis)
{
    vis[u]=1;
    split.push_back(u);
    if(node_in_cycle.find(u)!=node_in_cycle.end())
    {
        done=true;
        return;
    }
    for(auto x:tree[u])
    {
        if(!vis[x])
        {
            find_cycle_split(x,tree,node_in_cycle,split,done,vis);
            if(done)
                return;
        }
    }
    if(done)
        return;
    split.pop_back();
}

/** checks if a test point is inside the polygon **/
int pnpoly(vector<long double> vertx, vector<long double> verty, long double testx, long double testy)
{
  int nvert=(int)vertx.size(),i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

/** adds all points in the spanning tree which are inside the polygon formed by the cycle, to the cycle**/

/** the main function that identifies all the blocks **/

/**
void process_roundabouts()
{
    long long cnt_round=0;

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
    //cout<<cnt_round<<" roundabouts\n";
    //way_nodes.size()<<" "<<way_tags.size()<<"\n";
    cout<<cnt_round<<" ";
}
**/

/** ref: https://stackoverflow.com/questions/4681737/how-to-calculate-the-area-of-a-polygon-on-the-earths-surface-using-python/4683144 **/
long double block_area(vector<long long> nodes)
{
    long double pi = acos(-1.0), rad_earth = 6378100.0;
    long double lat_dist = rad_earth*pi/180.0;
    vector<long double> X,Y;    //X,Y coordinates
    for(auto v:nodes)
    {
        long double lat=(find_lat(v));
        long double lon=(find_lon(v));
        X.push_back(lat*lat_dist);
        Y.push_back(lon*lat_dist*cos(lat*pi/180.0));
    }
    long double area=0.0;
    long long N=X.size();
    for(int i=0;i<N;i++)
    {
        area+=(X[i]*Y[(i+1)%N]-Y[i]*X[(i+1)%N]);
    }
    area/=2.0;
    area=abs(area);
    return area;
}

/** Given the list of undirected edges (given in both directions) , this function finds all the regions created by these edges **/
void find_blocks()
{
    long double pi = acos(-1.0);
    vector<pair<pair<long long,long double>,long long>> list_with_angle;
    sort(list_of_edges.begin(),list_of_edges.end());
    //cout<<"Edges:\n";
    for(auto edge:list_of_edges)
    {
        long long U=edge.first,V=edge.second;
        //cout<<U<<" "<<V<<':';
        long double X1 = find_lon(U),Y1=find_lat(U),X2=find_lon(V),Y2=find_lat(V);
        long double angle = atan2(Y2-Y1,X2-X1);
        if(angle<0.0)
        {
            angle+=(2*pi);
        }
        //cout<<angle<<'\n';
        list_with_angle.push_back(make_pair(make_pair(U,angle),V));
    }
    sort(list_with_angle.begin(),list_with_angle.end());

    vector<pair<pair<long long,long long>,long long>> wedge_list;

    int i=0,j=0;
    while(i<list_with_angle.size())
    {
        long long cur = list_with_angle[i].first.first;
        j=i;
        //cout<<cur<<'\n';
        vector<pair<pair<long long,long double>,long long>> cur_list;
        while((j<list_with_angle.size()) and (list_with_angle[j].first.first==cur))
        {
            cur_list.push_back(list_with_angle[j]);
            j++;
        }
        for(int k=0;k<cur_list.size()-1;k++)
        {
            wedge_list.push_back(make_pair(make_pair(cur_list[k+1].second,cur),cur_list[k].second));
        }
        if(cur_list.size()>=2)
        {
            wedge_list.push_back(make_pair(make_pair(cur_list[0].second,cur),cur_list.back().second));
        }
        i=j;
    }
    sort(wedge_list.begin(),wedge_list.end());
    bool used[(long long)wedge_list.size()];
    memset(used,false,sizeof(used));
    vector<vector<long long>> region_list;
    for(int i=0;i<wedge_list.size();i++)
    {
        //cout<<wedge_list[i].first.first<<" "<<wedge_list[i].first.second<<" "<<wedge_list[i].second<<'\n';
        if(!used[i])
        {
            used[i]=true;
            long long v1=wedge_list[i].first.first,v2=wedge_list[i].first.second,v3=wedge_list[i].second;
            long long cur1=v2,cur2=v3;
            vector<long long> region;
            region.push_back(v3);
            bool found=true;
            while(true)
            {
                pair<pair<long long,long long>,long long> value = make_pair(make_pair(cur1,cur2),0LL);
                auto it1 = lower_bound(wedge_list.begin()+i,wedge_list.end(),value,[](const pair<pair<long long,long long>,long long> & a,const pair<pair<long long,long long>,long long> & b){
                return a.first.first<b.first.first;});
                auto it2 = upper_bound(wedge_list.begin()+i,wedge_list.end(),value,[](const pair<pair<long long,long long>,long long> & a,const pair<pair<long long,long long>,long long> & b){
                return a.first.first<b.first.first;});
                if((it1==wedge_list.end()) or ((*it1).first.first!=value.first.first))
                {
                    found=false;
                    break;
                }
                auto it = lower_bound(it1,it2,value,[](const pair<pair<long long,long long>,long long> & a,const pair<pair<long long,long long>,long long> & b){
                return a.first.second<b.first.second;});
                if((it==it2) or ((*it).first.second!=value.first.second))
                {
                    found=false;
                    break;
                }
                cur1=cur2;
                cur2=(*it).second;
                region.push_back(cur2);
                if(cur1==v1 and cur2==v2)
                {
                    break;
                }
            }
            if(found)
            {
                region_list.push_back(region);
            }
        }
    }
    long double avg_block_area=0.0;
    //cout<<"\n";
    //cout<<region_list.size()<<" Blocks found:\n";
    cout<<region_list.size()<<" ";
    int iter=0;
    for(auto region:region_list)
    {
        //cout<<"\nBlock "<<++iter<<":\n";
        //cout<<"No. of nodes:"<<region.size()<<'\n';
        //cout<<"Nodes:";

        /*for(auto u:region)
        {
            cout<<u<<" ";
        }
        cout<<'\n';*/
        long double area=block_area(region);
        //cout<<"Area:"<<area<<" sq meters";
        //cout<<'\n';
        avg_block_area+=area;
    }
    avg_block_area/=(long double)region_list.size();
    //cout<<"\nAverage Block Area:"<<avg_block_area<<" sq meters\n";
    cout<<avg_block_area<<" ";
}


// long long temp_count=0;
void find_intersections()
{
    long long count_of_intersections[max_intersections]={0};
    long long count3=0,count4=0,count_others=0;
    // cerr<<"intersections start";
    // cout<<node_store.size()<<endl;
    for(int i=0;i<node_store.size();i++)
    {
        /**
        if(!insignificant_node[i])
        {
            count_of_intersections[degree[i]]++;
            if(degree[i]==3)count3++;
            else if(degree[i]==4)count4++;
            else if(degree[i]>4)count_others++;
        }
        else
        {
            temp_count++;
            //cout<<node_data[i]["id"]<<'\n';

        }
        **/
        // cout<<i<<" "<<degree[i]<<endl;
        count_of_intersections[degree[i]]++;
        if(degree[i]==3)count3++;
        else if(degree[i]==4)count4++;
        else if(degree[i]>4)count_others++;
    }
    // cout<<temp_count<<'\n';
    /*
    cout<<"Total Nodes:"<<node_data.size()<<'\n';

    for(int i=0;i<max_intersections;i++)
    {
        if(count_of_intersections[i])
        cout<<"Number of "<<i<<"-way intersections is "<<count_of_intersections[i]<<endl;
    }*/
    // cerr<<"intersections done"<<endl;

    cout<<count3<<","<<count4<<","<<count_others<<",";
}

void walkability_ratio(long long n)
{
    // cerr<<"walk start!"<<endl;
    long double avg_ratio=0.0;
    long long N = node_data[n].size();
    if(N>walk_ratio_index)
    {
        for(int i=0;i<epochs;i++)
        {
        mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
        vector<pair<long long,long long>> node_list;


        map<pair<long long,long long>,long long> used;
        used.clear();
        while(node_list.size()<walk_ratio_index)
        {
            long long X=uniform_int_distribution<long long>(0, N-1)(rng);
            long long Y=uniform_int_distribution<long long>(0, N-1)(rng);
            string ID_string1 =node_data[n][X]["id"];
            ID_string1=ID_string1.substr(1,ID_string1.size()-2);
            long long ID1 = stoll(ID_string1);
            string ID_string2 =node_data[n][Y]["id"];
            ID_string2=ID_string2.substr(1,ID_string2.size()-2);
            long long ID2 = stoll(ID_string2);
            X = node_id_map[ID1];
            Y = node_id_map[ID2];
            // cout<<X<<" "<<Y<<endl;
            if(X==Y)continue;
            pair<long long,long long> cur_pair=make_pair(X,Y);
            if(used.find(cur_pair)!=used.end()) continue;
            used[make_pair(X,Y)]=1;
            /**
            if(!insignificant_node[X] and !insignificant_node[Y])
            {
                node_list.emplace_back(X,Y);
            }
            **/
            node_list.emplace_back(X,Y);
        }
        //cout<<"Epoch No "<<i<<endl;
        /**
        cout<<"Nodes used for calculating walkability ratio:\n";
        for(auto u:node_list)
        {
            cout<<node_data[u]["id"]<<" ";
        }
        cout<<endl;
        **/
        // cout<<"Points selected!"<<endl;
        long double ratio = 0.0;
        for(int i=0;i<walk_ratio_index;i++)
        {
            vector<long double> dist;
            vector<long long> pre;
            dijkstra(n,current_map[node_list[i].first],dist,pre,current_map[node_list[i].second]);
            long double beeline_dist = find_distance(find_lat(node_list[i].first),find_lon(node_list[i].first),find_lat(node_list[i].second),find_lon(node_list[i].second));
            //cout<<dist[node_list[i].second]<<endl;
            ratio = ratio + beeline_dist/(dist[current_map[node_list[i].second]]);
            if(!dist[current_map[node_list[i].second]]){
                cerr<<"issue with dijkstra"<<endl;
                break;
            }
        }
        // cout<<"done"<<endl;
        ratio/=walk_ratio_index;
        // long double dist[walk_ratio_index+1][walk_ratio_index+1];
        // for(int i=0;i<walk_ratio_index;i++)
        // {
        //  for(int j=0;j<walk_ratio_index;j++)
        //  {
        //      if(edge_list.find(make_pair(node_list[i],node_list[j]))!=edge_list.end())
        //          dist[i][j]=edge_list[make_pair(node_list[i],node_list[j])];
        //      else
        //          dist[i][j]=INF;
        //  }
        // }
        // for(int i=0;i<walk_ratio_index;i++)
        // {
        //  dist[i][i]=0;
        // }
        // floyd_warshall(dist);
        // long double ratio=0.0,pair_count=0.0;
        // for(int i=0;i<walk_ratio_index;i++)
        // {
        //  for(int j=0;j<walk_ratio_index;j++)
        //  {
        //      if(i!=j)
        //      {

     //                if(node_list[i]!=node_list[j])
     //                {
     //                    ratio = ratio + (beeline_dist)/dist[i][j];        //do this only when different nodes are encountered
     //                    pair_count=pair_count+1.0;
     //                }

        //      }
        //  }
        // }
        // ratio/=(pair_count);
            avg_ratio+=ratio;
        }
    avg_ratio/=epochs;
    //cout<<"\nWalkability Ratio:\n";
    cout<<avg_ratio<<",";
    }
    else
    {
        cout<<"NA,";
    }
    // cerr<<"walk done";
    // cerr<<endl;
}


void walkability_ratio_v2(long long n){
    // cerr<<"walk_v2 enter "<<n<<endl;
    long double avg_ratio = 0.0;
    long long N = node_data[n].size();
    for(int i=0;i<epochs;i++){
        mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
        vector<pair<long long,long long>> node_list;
        long long max_scale = 1e9;
        int iter=0;
        long double ratio = 0.0;
        while(iter<walk_ratio_index){
            // cerr<<"in loop"<<endl;
            long double X1=uniform_int_distribution<long long>(0, max_scale-1)(rng);
            long double Y1=uniform_int_distribution<long long>(0, max_scale-1)(rng);
            X1/=max_scale;
            Y1/=max_scale;
            X1 = minlat + (n%rows)*0.01 + 0.01*X1;
            Y1 = minlon + (n/rows)*0.01 + 0.01*Y1;
            long double X2=uniform_int_distribution<long long>(0, max_scale-1)(rng);
            long double Y2=uniform_int_distribution<long long>(0, max_scale-1)(rng);
            X2/=max_scale;
            Y2/=max_scale;
            X2 = minlat + (n%rows)*0.01 + 0.01*X2;
            Y2 = minlon + (n/rows)*0.01 + 0.01*Y2;
            long double beeline_dist = find_distance(X1,Y1,X2,Y2), path_dist=0.0;
            long double min_dist = 1e18;
            long long min_A = -1;
            // cerr<<"second loop "<<node_data[n].size()<<" " <<endl;
            for(int j=0;j<node_data[n].size();j++){
                // cerr<<"intheloop"<<endl;
                string ID_string1 =node_data[n][j]["id"];
                ID_string1=ID_string1.substr(1,ID_string1.size()-2);
                long long ID1 = stoll(ID_string1);
                long long A = node_id_map[ID1];
                long double cur_dist = find_distance(X1,Y1,find_lat(A),find_lon(A));
                if(cur_dist<min_dist){
                    min_dist = cur_dist;
                    min_A = A;
                }
            }
            // cerr<<"here dumbo"<<endl;
            if(min_A==-1){
                // cerr<<"nothing found"<<endl;
                break;
            }
            path_dist+=min_dist;
            long long first_point,second_point;
            first_point = min_A;
            min_dist = 1e18;
            min_A = -1;
            for(int j=0;j<node_data[n].size();j++){
                // cerr<<"here haha"<<endl;
                string ID_string1 =node_data[n][j]["id"];
                ID_string1=ID_string1.substr(1,ID_string1.size()-2);
                long long ID1 = stoll(ID_string1);
                long long A = node_id_map[ID1];
                long double cur_dist = find_distance(X2,Y2,find_lat(A),find_lon(A));
                if(cur_dist<min_dist){
                    min_dist = cur_dist;
                    min_A = A;
                }
            }
            second_point = min_A;
            path_dist+=min_dist;
            // cerr<<"path_dist:"<<path_dist<<endl;
            vector<long double> dist;
            vector<long long> pre;
            // cerr<<"before dijkstra"<<endl;
            // cerr<<current_map[first_point]<<" "<<current_map[second_point]<<" "<<current_map.size()<<" " <<endl;

            dijkstra(current_map.size()+1,current_map[first_point],dist,pre,current_map[second_point]);
            // cerr<<"path_dist:"<<path_dist<<endl;
            // cerr<<"after dijkstra"<<endl;
            path_dist+=dist[current_map[second_point]];
            // cerr<<"path_dist:"<<path_dist<<endl;
            ratio = ratio + (beeline_dist/path_dist);
            iter++;
            // cerr<<"beeline_dist:"<<beeline_dist<<endl;
            // cerr<<"path_dist:"<<path_dist<<endl;
        }
        ratio/=walk_ratio_index;
        avg_ratio+=ratio;

    }
    avg_ratio/=epochs;
    cout<<avg_ratio<<'\n';
    // cerr<<"walk_v2 out "<<n<<endl;
}

/** note that lon acts like the x-coordinate while
lat acts like the y-coordinate**/
set<pair<long long, long long>> set_of_all_edges;
void total_road_length(long long n)
{
    long double road_length=0.0;
    for(auto u:edge_list_all[n])
    {
        long long U,V;
        tie(U,V) = u.first;
        string lat_val1 = node_store[U]["lat"],lat_val2 = node_store[V]["lat"];
        long double X1=stold(lat_val1.substr(1,lat_val1.size()-2)),X2=stold(lat_val2.substr(1,lat_val2.size()-2));
        string lon_val1 = node_store[U]["lon"],lon_val2 = node_store[V]["lon"];
        long double Y1=stold(lon_val1.substr(1,lon_val1.size()-2)),Y2=stold(lon_val2.substr(1,lon_val2.size()-2));
        long long j1 = floor(100.0*(X1-minlat)), i1 = floor(100.0*(Y1-minlon));
        long long j2 = floor(100.0*(X2-minlat)), i2 = floor(100.0*(Y2-minlon));
        long long n1 = i1*rows + j1, n2 = i2*rows + j2;
        set_of_all_edges.insert(make_pair(U,V));
        if((n1==n) and (n2==n))
            road_length+=(find_distance(find_lat(u.first.first),find_lon(u.first.first),find_lat(u.first.second),find_lon(u.first.second)))/2.0;
        else{
            long double Xn,Yn;
            if(n1==n){
                n1 = n2;
                Xn = Y1;
                Yn = X1;
                swap(X2,Y2);
            }
            else{
                Xn = Y2;
                Yn = X2;
                X2 = Y1;
                Y2 = X1;
            }
            long double eps = 0.000000000001;
            long double min_X = floor(100.0*Xn)/100.0, max_X = ceil(100*(Xn+eps))/100.0, min_Y = floor(100.0*Yn)/100.0, max_Y = ceil(100.0*(Yn+eps))/100.0;
            if(n1==(n-1)){
                long double X,Y=min_Y;
                X = ((X2-Xn)*(Y-Yn))/(Y2-Yn) + Xn;
                road_length+=find_distance(Yn,Xn,Y,X)/2;
            }
            else if(n1==(n+1)){
                long double X,Y=max_Y;
                X = ((X2-Xn)*(Y-Yn))/(Y2-Yn) + Xn;
                road_length+=find_distance(Yn,Xn,Y,X)/2;   
                if(n==13){
                    cerr<<"Interpolate:\n";
                    cerr<<X<<" "<<Y<<'\n';
                }
            }
            else if(n1==(n-rows)){
                long double X = min_X,Y;
                Y = ((Y2-Yn)*(X-Xn))/(X2-Xn) + Yn;
                road_length+=find_distance(Yn,Xn,Y,X)/2;   
            }
            else if(n1==(n+rows)){
                long double X = max_X,Y;
                Y = ((Y2-Yn)*(X-Xn))/(X2-Xn) + Yn;
                road_length+=find_distance(Yn,Xn,Y,X)/2;   
            }

        }
    }
    // cout<<"\nTotal Road Length:\n";
    cout<<road_length<<",";
    // cerr<<"road done"<<endl;
}
/**
This function will first determine the grid no. for which the graph is being
created using the minlat,maxlat,etc. values.
It creates the
**/
int revisions[5000];
int find_year(string &s){
    int y=0;
    for(auto u:s){
        y*=10;
        y+=(u-'0');
    }
    return y;
}
void add_edges(long long n, set<pair<long long,long long>> &set_of_edges, long long &iterator)
{
    /**
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
                string lat_val1 = node_data[U]["lat"],lat_val2 = node_data[V]["lat"];
                long double X1=stold(lat_val1.substr(1,lat_val1.size()-2)),X2=stold(lat_val2.substr(1,lat_val2.size()-2));
                string lon_val1 = node_data[U]["lon"],lon_val2 = node_data[V]["lon"];
                long double Y1=stold(lon_val1.substr(1,lon_val1.size()-2)),Y2=stold(lon_val2.substr(1,lon_val2.size()-2));
                long double dist = find_distance(X1,Y1,X2,Y2);
                edge_list[E1]=dist;
                edge_list[E2]=dist;

                adj_dist[U].emplace_back(V,dist);
                adj_dist[V].emplace_back(U,dist);

                list_of_edges.push_back(E1);
                list_of_edges.push_back(E2);
            }
        }
    }
    **/


    /** add mapped vertices to make it fast **/


    // cerr<<"add_edge "<<n<<endl;
    long long it = iterator;
    // if(n==3)
    //     cerr<<n<<'\n';
    for(auto u:node_data[n]){
        string ID_string =u["id"];
        ID_string=ID_string.substr(1,ID_string.size()-2);
        long long ID = stoll(ID_string);
        // if(n==3){
        //     cerr<<ID<<'\n';
        // }
        // cerr<<u["timestamp"]<<'\n';
        // string year, timestamp = u["timestamp"];
        // year = timestamp.substr(1,4);
        // revisions[find_year(year)]++;
        ID = node_id_map[ID];
        current_map[ID] = it++;
    }
    iterator=it;

    for(auto E:edge_list_all[n]){
        long long U = E.first.first, V = E.first.second;
        if(set_of_edges.find(E.first)!=set_of_edges.end())
            continue;
        set_of_edges.insert(E.first);
        // cerr<<U<<" "<<V<<endl;
        // degree[U]++;
        if((current_map.find(U)==current_map.end()) and (current_map.find(V)==current_map.end()))
        {
            cerr<<"error!!!"<<endl;
            break;
        }
        if(current_map.find(V)!=current_map.end()){
            V = current_map[V];

            degree[V]++;
        

            if(current_map.find(U)!=current_map.end()){
                U = current_map[U];
                adjacency_list[U].push_back(V);

                adj_dist[U].emplace_back(V,E.second);

            }
            

            // adjacency_list[V].push_back(U);
        }
        // adj_dist[V].emplace_back(U,E.second);
    }
    // cerr<<"added edges"<<endl;

}


void refresh_graph(long long n){
    // cerr<<"refresh"<<endl;
    for(auto u:node_data[n]){
        string ID_string =u["id"];
        ID_string=ID_string.substr(1,ID_string.size()-2);
        long long ID = stoll(ID_string);
        ID = node_id_map[ID];
        ID = current_map[ID];
        adjacency_list[ID].clear();
        adj_dist[ID].clear();
        degree[ID]=0;
    }
    // current_map.clear();
    // cerr<<"refresh_graph done"<<endl;
}

void way(ifstream & inp_handle)
{
    vector<string> way_strings = split_string(s,' ');

    
    // No Need of way_data now!


    /*
    map<string,string> way_data_map;
    way_data_map.clear();
    for(int i=1;i<way_strings.size()-1;i++)                        //first and last elements are not data elements
    {
        vector<string> temp = split_string(way_strings[i],'=');
        if(temp.size()>=2)
        way_data_map[temp[0]]=temp[1];
    }

    string ID_string =way_data_map["id"];
    ID_string=ID_string.substr(1,ID_string.size()-2);
    long long ID = stoll(ID_string);

    way_id_map[ID]=(long long)way_data.size();                   //for almost O(1) access
    way_data.push_back(way_data_map);*/

    // store node data for every way
    vector<map<string,string>> cur_way_node_data,cur_way_tag_data;
    // cout<<s<<endl;
    while(getline(inp_handle,s))
    {
        // cout<<s<<endl;
        s=ignore_leading_space(s);
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
                cur_way_tag_data_map[temp[0]]=temp[1];
            }
            cur_way_tag_data.push_back(cur_way_tag_data_map);
        }
        else
            break;

    }


    /** No need of way_nodes or way_tags in the graph.
        We can directly store the edges from here itself!
    **/

    for(int i=0;i<cur_way_node_data.size()-1;i++)
    {
        string ID_string1 = cur_way_node_data[i]["ref"];
        string ID_string2 = cur_way_node_data[i+1]["ref"];
        // cerr<<"Edge:\n"<<ID_string1<<" "<<ID_string2<<'\n';
        ID_string1=ID_string1.substr(1,ID_string1.size()-2);
        ID_string2=ID_string2.substr(1,ID_string2.size()-2);
        long long ID1 = stoll(ID_string1), ID2= stoll(ID_string2);
        if((node_id_map.find(ID1)==node_id_map.end()) or (node_id_map.find(ID2)==node_id_map.end())) {
            continue;
        }
        long long U=node_id_map[ID1],V=node_id_map[ID2];
        string lat_val1 = node_store[U]["lat"],lat_val2 = node_store[V]["lat"];
        long double X1=stold(lat_val1.substr(1,lat_val1.size()-2)),X2=stold(lat_val2.substr(1,lat_val2.size()-2));
        string lon_val1 = node_store[U]["lon"],lon_val2 = node_store[V]["lon"];
        long double Y1=stold(lon_val1.substr(1,lon_val1.size()-2)),Y2=stold(lon_val2.substr(1,lon_val2.size()-2));
        // cerr<<"Edge:\n";
        // cerr<<U+1<<" "<<V+1<<'\n';
        if((X1<minlat) or (X1>maxlat) or (Y1<minlon) or (Y1>maxlon)){
            // cout<<"removed"<<endl;
            continue;
        }
        if((X2<minlat) or (X2>maxlat) or (Y2<minlon) or (Y2>maxlon)){
            // cout<<"removed"<<endl;
            continue;
        }
        long long j1 = floor(100.0*(X1-minlat)), i1 = floor(100.0*(Y1-minlon));
        long long j2 = floor(100.0*(X2-minlat)), i2 = floor(100.0*(Y2-minlon));
        long long n1 = i1*rows + j1, n2 = i2*rows + j2;

        pair<long long,long long> E1=make_pair(U,V),E2=make_pair(V,U);

        long double dist = find_distance(X1,Y1,X2,Y2);
        edge_list_all[n1][E1]=dist;
        edge_list_all[n1][E2]=dist;
        if(n1==3){
            // cerr<<"Edge test:\n";
            // cerr<<U+1<<" "<<V+1<<'\n';
            // cerr<<n2<<'\n';
        }
        if(n1!=n2){
            edge_list_all[n2][E1]=dist;
            edge_list_all[n2][E2]=dist;
        }
        // if((edge_list_all[n1].find(E1)==edge_list_all[n1].end()) and (edge_list_all[n1].find(E2)==edge_list_all[n1].end()) and (U!=V))
        // {

        //     adjacency_list[U].push_back(V);
        //     adjacency_list[V].push_back(U);
        //     degree[U]++;
        //     degree[V]++;
        //     long double dist = find_distance(X1,Y1,X2,Y2);
        //     edge_list_all[n1][E1]=dist;
        //     edge_list_all[n1][E2]=dist;

        //     adj_dist[U].emplace_back(V,dist);
        //     adj_dist[V].emplace_back(U,dist);

        //     list_of_edges.push_back(E1);
        //     list_of_edges.push_back(E2);
        // }
    }
    // way_nodes.push_back(cur_way_node_data);
    // way_tags.push_back(cur_way_tag_data);
}
void node()
{

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
    // cout<<ID<<endl;

    string lat_val = node_data_map["lat"],lon_val = node_data_map["lon"];
    long double X=stold(lat_val.substr(1,lat_val.size()-2)),Y=stold(lon_val.substr(1,lon_val.size()-2));
    if((X<minlat) or (X>maxlat) or (Y<minlon) or (Y>maxlon)){
        // cout<<"removed"<<endl;
        return;
    }
    long long j = floor(100.0*(X-minlat)), i = floor(100.0*(Y-minlon));
    long long n = i*rows + j;
    node_data[n].push_back(node_data_map);
    // cerr<<"Node:\n"<<ID<<" "<<(long long)node_store.size()<<'\n';
    node_id_map[ID]=(long long)node_store.size();                 //for almost O(1) access
    node_store.push_back(node_data_map);
}

/**
Some useful conversions from string
stold -> string to long double
stoll -> string to long long
**/

vector<vector<int>> find_urban_extent(vector<vector<int>> bu_matrix)
{
    vector<vector<int>> label = bu_matrix;  //label to indicate the type of the current pixel
    int N = bu_matrix.size(), M = bu_matrix[0].size();
    //ri -> 25 rj -> 20
    int ri = 25, rj = 20;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            int built = 0, tot = 0;
            for(int u=i-ri ; u<=i+ri;u++)
            {
                for(int v = j-rj; v<=j+rj; v++)
                {
                    if(u<0 or u>=N or v<0 or v>=M)
                        continue;
                    if((400.0*((double)u*u) + 625.0*((double)v*v) - 250000.0)<=0)
                    {
                        built+=(bu_matrix[u][v]==1);
                    }
                    tot++;
                }
            }
            long double perc = ((long double)(built))/tot;
            if(perc>=0.5)
                label[i][j] = 0;    //urban
            else if(perc>=0.25)
                label[i][j] = 1;    //suburban
            else
                label[i][j] = 2;    //rural
        }
    }
    //labelled all pixels
}
long double latlon[641][4] = { {18.68,77.76,19.92,79.98} ,{26.74,77.43,27.41,78.85},{21.97,71.61,23.51,72.84},{18.33,73.62,19.99,75.59},{23.27,92.62,24.41,93.21},{25.65,73.91,26.98,75.36},{20.28,76.69,21.25,77.62},{9.08,76.28,9.89,76.69},{27.58,77.48,28.17,78.62},{21.92,74.03,22.65,74.75},{24.81,81.52,25.75,82.35},{29.43,79.04,29.98,80.07},{27.06,76.12,28.22,77.21},{30.07,76.54,30.58,77.21},{26.16,82.22,26.65,83.13},{20.54,76.63,21.77,78.45},{20.69,70.64,22.04,72.1},{31.48,74.49,32.05,75.4},{22.15,72.33,22.74,73.23},{13.68,76.76,15.23,78.47},{33.35,75.06,34.23,75.71},{27.62,96.28,28.46,97.42},{20.53,84.27,21.68,85.4},{22.65,81.12,23.42,82.2},{25.94,87.04,26.59,87.7},{10.88,78.94,11.42,79.51},{24.23,77.48,25.01,78.28},{26.36,79.23,26.95,79.76},{24.48,84.0,25.13,84.75},{19.38,74.59,20.67,75.89},{25.64,82.67,26.45,83.51},{33.66,74.39,34.07,74.93},{15.81,74.99,16.77,76.33},{29.68,79.47,30.32,80.16},{28.78,77.13,29.3,77.51},{27.06,81.04,28.4,81.77},{26.46,90.82,26.83,91.79},{21.32,79.51,22.4,81.05},{20.15,82.69,21.07,83.68},{20.77,86.36,21.98,87.49},{25.56,83.68,26.19,84.64},{27.05,82.02,27.84,82.77},{23.83,71.26,24.71,73.03},{24.98,80.1,25.91,81.04},{34.31,74.37,34.87,75.43},{12.35,77.18,13.5,77.97},{12.66,77.33,13.23,77.84},{24.54,86.49,25.12,87.17},{22.63,86.61,23.64,87.77},{23.06,73.96,23.92,74.77},{26.52,80.93,27.36,81.73},{20.73,82.64,21.74,83.91},{33.86,74.02,34.4,74.82},{24.4,76.22,25.43,77.42},{22.93,86.8,23.88,88.42},{28.02,79.0,28.89,79.76},{24.65,70.09,26.51,72.85},{30.14,75.26,30.64,75.73},{26.09,90.66,26.66,91.31},{21.37,74.43,22.13,75.63},{18.69,81.29,20.19,82.25},{26.54,82.22,27.13,82.98},{29.78,74.63,30.58,75.38},{20.39,83.57,20.9,84.8},{25.25,85.74,25.78,86.52},{15.36,74.09,16.96,75.46},{14.56,75.66,15.83,77.17},{21.36,77.01,22.4,78.55},{20.74,86.27,21.24,86.98},{25.05,86.64,25.51,87.55},{20.64,79.45,21.6,80.09},{26.71,76.88,27.82,77.76},{21.42,72.5,22.25,73.5},{21.0,71.39,22.35,72.37},{25.01,74.01,25.96,75.46},{25.9,78.21,26.79,79.14},{28.38,75.48,29.08,76.47},{25.16,84.28,25.74,84.85},{23.09,77.17,23.9,77.65},{18.45,74.81,19.44,76.74},{17.58,76.69,18.48,77.66},{18.14,80.25,19.4,81.23},{16.15,75.33,17.48,76.47},{29.03,77.99,29.79,78.94},{27.19,71.9,29.05,74.36},{21.71,81.48,23.12,82.48},{31.21,76.39,31.6,76.93},{23.54,87.09,24.58,88.03},{24.29,93.71,24.74,93.89},{23.42,85.58,23.97,86.45},{26.18,90.38,26.53,90.86},{27.65,78.28,28.48,79.5},{28.07,77.64,28.71,78.49},{19.84,75.93,21.29,76.83},{25.0,75.26,25.88,76.33},{21.07,75.95,21.6,76.79},{25.26,83.77,25.75,84.4},{24.37,92.42,25.14,93.26},{28.62,77.18,28.66,77.26},{32.18,75.8,33.21,76.88},{29.93,79.08,31.07,80.1},{29.05,79.8,29.53,80.32},{23.0,93.01,24.08,93.44},{11.58,76.4,12.31,77.78},{24.71,83.01,25.53,83.56},{23.84,93.74,24.64,94.41},{30.67,76.7,30.8,76.85},{19.46,78.8,20.73,80.0},{26.88,95.63,27.66,97.17},{23.68,84.44,24.53,85.35},{12.99,80.2,13.14,80.31},{24.1,78.99,25.42,80.43},{21.46,78.25,22.82,79.41},{13.22,77.36,13.96,78.21},{12.91,75.08,13.9,76.36},{26.4,90.35,26.9,90.96},{13.57,76.03,15.04,77.03},{24.89,80.69,25.54,81.56},{24.22,74.11,25.22,75.82},{12.62,78.07,13.99,80.05},{23.94,92.97,24.61,93.88},{27.41,73.86,29.0,75.67},{10.22,76.66,11.41,77.3},{11.15,78.88,11.9,79.81},{13.72,77.95,15.24,79.47},{19.96,84.85,20.7,86.43},{20.05,72.91,20.36,73.22},{17.79,80.92,19.34,81.97},{25.17,88.19,25.59,89.01},{12.46,74.67,13.67,75.67},{20.37,72.8,20.47,72.9},{23.15,79.06,24.45,79.96},{25.87,85.67,26.45,86.41},{26.45,87.99,27.22,88.89},{26.2,91.75,26.68,92.39},{32.76,72.53,37.08,77.49},{25.55,78.22,26.29,78.84},{13.81,75.4,14.93,76.53},{21.13,84.48,21.74,85.22},{29.94,77.58,30.98,78.31},{24.04,86.46,24.62,87.07},{26.09,83.48,26.76,84.19},{22.29,75.9,23.32,77.13},{23.41,91.75,24.24,92.17},{20.05,81.41,21.02,82.17},{23.62,86.11,24.06,86.83},{22.02,74.47,23.14,75.71},{11.77,77.67,12.51,78.75},{15.04,74.73,15.7,75.56},{26.36,77.23,26.95,78.27},{27.21,94.21,27.88,95.52},{20.49,84.87,21.19,86.03},{25.48,89.7,26.39,90.48},{20.63,73.85,21.63,75.21},{28.47,95.25,29.46,96.62},{27.09,94.57,27.71,95.49},{10.02,77.27,10.83,78.34},{22.44,80.48,23.37,81.74},{20.7,70.87,20.77,71.16},{32.85,75.29,33.4,76.11},{22.5,73.8,23.33,74.48},{23.98,86.89,24.64,87.7},{23.33,73.36,24.02,74.38},{20.38,80.81,22.02,81.94},{25.65,93.54,25.96,94.0},{28.58,77.25,28.67,77.35},{27.14,88.44,27.42,88.93},{25.41,90.12,26.02,91.03},{16.31,81.51,18.01,82.6},{26.91,92.42,27.98,93.39},{25.12,91.36,25.65,92.14},{21.54,76.02,22.42,77.22},{27.74,94.75,28.57,95.57},{9.79,76.17,10.3,76.84},{11.02,76.84,11.95,77.93},{27.31,78.18,27.78,79.28},{26.43,78.74,27.01,79.34},{26.41,81.54,26.92,82.48},{28.16,77.15,28.51,77.55},{30.37,74.47,30.83,75.06},{27.16,79.12,27.71,79.75},{29.25,75.22,29.82,75.96},{30.41,76.08,30.88,76.59},{25.44,80.22,26.23,81.34},{26.88,78.2,27.51,78.83},{29.95,73.89,31.16,75.11},{14.94,75.27,15.89,76.05},{18.75,83.8,19.64,84.45},{34.17,74.75,34.46,75.63},{22.99,72.34,23.57,73.03},{28.71,72.65,30.2,74.31},{18.98,84.13,20.28,85.18},{18.68,79.74,20.83,80.9},{23.56,83.33,24.54,84.07},{29.45,78.19,30.31,79.23},{28.08,77.3,28.68,77.76},{24.26,84.29,25.07,85.4},{28.55,77.2,28.93,78.21},{25.31,83.06,25.89,83.96},{23.89,85.67,24.78,86.58},{25.88,90.12,26.24,91.1},{24.5,87.05,25.23,87.52},{25.81,93.29,26.84,94.18},{26.78,81.51,27.44,82.62},{20.65,79.8,21.63,80.68},{26.2,83.91,26.64,84.91},{26.22,83.07,27.14,83.67},{16.71,76.07,17.77,77.69},{22.71,84.04,23.61,85.02},{23.89,76.81,25.11,77.74},{15.74,79.2,16.83,80.91},{31.59,74.89,32.58,75.95},{28.2,76.66,28.54,77.24},{25.56,77.67,26.35,78.9},{24.13,92.42,24.88,92.78},{31.42,76.3,31.9,76.73},{25.4,79.36,26.16,80.35},{28.78,73.8,29.96,75.53},{22.22,87.85,22.78,88.37},{21.91,76.78,22.58,77.5},{26.89,79.69,27.78,80.83},{29.54,77.71,30.25,78.34},{12.51,75.56,13.55,76.64},{27.27,77.88,27.83,78.53},{14.28,75.02,15.16,75.82},{23.66,85.02,24.36,85.93},{19.07,76.68,20.02,77.49},{28.9,75.26,29.58,76.31},{22.22,77.21,22.99,78.7},{31.13,75.47,32.08,76.35},{22.6,87.51,23.23,88.51},{17.3,78.39,17.48,78.54},{9.27,76.63,10.36,77.41},{24.55,93.77,25.07,94.14},{24.55,93.07,24.86,93.25},{22.33,75.42,23.09,76.25},{22.83,79.35,23.61,80.58},{19.97,86.02,20.4,86.79},{25.03,91.99,25.75,92.8},{26.44,74.92,27.86,76.29},{26.02,69.48,28.04,72.34},{20.58,85.69,21.17,86.63},{30.97,75.07,31.62,75.95},{25.77,78.93,26.44,79.96},{20.27,74.76,21.41,76.4},{19.28,75.59,20.65,76.53},{24.61,71.19,25.81,73.1},{26.26,88.4,27.0,89.88},{32.47,74.45,33.03,75.16},{21.7,68.94,22.96,70.66},{23.8,86.47,24.17,87.3},{24.37,85.83,25.14,86.62},{21.68,82.31,22.26,83.32},{22.29,83.4,23.25,84.4},{25.39,82.12,26.2,83.09},{22.39,74.16,23.24,75.01},{28.36,76.28,28.86,76.96},{23.76,75.46,24.87,76.95},{25.11,78.3,25.95,79.42},{21.57,83.43,22.04,84.39},{27.64,75.03,28.52,76.1},{25.85,71.79,27.62,73.87},{26.35,93.82,27.18,94.61},{20.73,69.94,21.68,71.22},{28.4,78.06,29.08,78.65},{22.73,68.19,24.69,71.72},{24.54,83.32,25.41,83.9},{29.51,75.95,30.2,76.76},{19.18,82.53,20.45,83.79},{25.72,90.94,26.57,91.82},{25.86,91.45,26.25,92.18},{12.23,79.56,13.14,80.27},{19.57,83.5,20.69,84.59},{31.69,75.59,32.47,77.08},{19.69,80.4,20.56,81.82},{26.77,79.32,27.23,80.02},{8.08,77.1,8.58,77.59},{11.67,75.17,12.3,75.94},{26.09,79.5,26.85,80.18},{25.92,79.94,26.96,80.57},{27.55,78.47,28.05,79.21},{31.12,74.95,31.65,75.91},{25.53,92.15,26.61,93.9},{32.89,75.36,34.86,77.53},{24.25,92.21,24.92,92.59},{17.98,78.52,19.07,80.34},{29.4,76.45,29.99,77.22},{10.55,77.75,11.09,78.59},{12.06,74.86,12.8,75.43},{32.28,75.19,32.85,75.94},{25.2,87.21,25.88,88.08},{23.31,79.82,24.14,80.96},{25.26,81.15,25.81,81.82},{20.3,86.25,20.8,87.1},{21.01,85.19,22.16,86.38},{25.25,86.28,25.73,86.86},{16.77,79.79,18.63,81.8},{22.5,72.51,23.29,73.58},{27.67,80.03,28.69,81.31},{19.68,84.94,20.43,86.08},{22.57,84.94,23.28,85.62},{31.1,77.75,32.09,79.01},{25.59,94.58,26.04,95.05},{25.93,87.61,26.56,88.3},{33.0,75.47,34.22,76.8},{25.93,88.64,26.55,89.87},{11.93,75.37,12.83,76.19},{24.1,85.08,24.82,85.91},{25.52,93.92,26.02,94.36},{26.13,89.83,26.9,90.43},{12.76,77.84,13.6,78.59},{23.95,92.53,24.52,92.9},{15.74,73.69,17.18,74.7},{22.5,88.27,22.63,88.41},{8.76,76.48,9.17,77.27},{15.14,75.77,16.01,76.81},{18.24,82.09,19.24,83.42},{22.04,82.14,23.0,83.12},{22.94,81.58,23.92,82.76},{24.54,75.62,25.84,76.58},{9.4,76.37,9.86,76.99},{11.12,75.54,11.8,76.14},{15.71,80.0,17.15,81.56},{12.14,77.47,12.89,78.67},{33.45,74.51,33.83,75.16},{31.34,76.94,32.42,77.87},{34.29,73.91,34.87,74.6},{14.9,76.97,16.17,78.93},{29.86,76.42,30.25,77.12},{27.61,92.67,28.37,94.02},{26.55,83.53,27.3,84.42},{31.75,76.38,33.25,78.67},{26.82,93.7,27.52,94.58},{24.97,85.9,25.34,86.39},{24.18,78.17,25.22,78.99},{23.32,83.97,24.05,84.97},{17.87,76.2,18.84,77.3},{21.98,92.51,22.76,92.97},{32.34,76.38,36.0,80.33},{23.28,84.4,23.68,84.95},{27.54,95.76,28.37,96.7},{26.37,94.67,26.8,94.92},{27.88,95.32,28.8,96.43},{27.34,93.52,27.96,94.36},{26.5,80.56,27.16,81.22},{30.56,75.36,31.02,76.37},{22.5,92.36,23.4,93.17},{25.44,86.61,26.12,87.12},{26.04,85.75,26.67,86.72},{9.56,77.46,10.31,78.47},{26.88,83.12,27.48,83.94},{20.82,82.0,21.55,83.28},{15.83,77.24,17.24,79.24},{11.69,75.53,11.76,75.56},{27.8,75.9,28.47,76.37},{23.04,71.94,24.09,72.87},{25.03,79.28,25.64,80.15},{26.93,78.71,27.47,79.44},{10.69,75.83,11.53,76.55},{24.66,87.76,25.54,88.47},{17.81,81.39,18.74,82.45},{23.26,92.26,24.25,92.68},{31.23,76.62,32.07,77.39},{22.2,79.95,23.17,81.2},{23.76,74.88,24.76,75.93},{12.22,76.33,13.05,77.33},{29.54,75.17,30.21,75.78},{26.06,91.96,26.52,92.56},{27.23,77.28,27.96,77.95},{25.8,83.31,26.28,83.81},{21.27,85.66,22.57,87.18},{17.42,77.44,18.29,79.14},{28.74,77.43,29.26,78.14},{27.65,76.85,28.33,77.34},{24.6,82.08,25.27,83.18},{30.49,74.9,31.1,75.42},{26.19,94.29,26.76,94.76},{26.43,94.78,27.04,95.24},{28.33,78.41,29.28,78.98},{25.9,77.12,26.87,78.54},{29.9,74.25,30.67,74.82},{18.85,72.79,19.05,72.91},{18.98,72.78,19.27,72.98},{24.95,86.3,25.5,86.74},{23.72,87.82,24.86,88.74},{29.18,77.1,29.71,78.14},{25.9,84.88,26.39,85.75},{11.74,75.91,12.66,77.13},{19.15,81.84,20.11,82.86},{22.87,88.14,24.2,88.8},{25.71,92.4,26.7,93.32},{26.41,73.11,27.71,75.37},{20.58,78.25,21.72,79.66},{28.98,78.84,29.61,79.98},{24.97,85.17,25.46,85.93},{26.14,91.22,26.58,91.57},{16.37,78.61,17.81,80.08},{11.01,77.68,11.59,78.49},{18.26,76.93,19.93,78.37},{21.0,73.58,22.03,74.77},{19.23,80.66,20.09,81.52},{21.39,73.3,22.08,73.99},{22.62,78.44,23.25,79.64},{19.59,73.25,20.87,74.94},{20.58,72.7,21.07,73.5},{24.52,85.27,25.11,86.06},{19.9,84.49,20.58,85.46},{24.24,74.71,25.05,75.61},{13.47,79.07,15.12,80.27},{28.58,77.17,28.64,77.26},{18.08,77.52,19.0,78.68},{28.65,77.17,28.79,77.25},{27.38,88.12,28.13,88.89},{21.55,88.33,23.25,89.1},{24.97,92.53,25.82,93.48},{28.66,77.21,28.78,77.34},{15.27,73.68,15.8,74.28},{23.65,91.91,24.53,92.34},{28.67,76.95,28.88,77.23},{19.99,82.33,21.09,82.86},{17.64,75.29,18.7,76.79},{24.24,87.39,24.83,87.92},{10.33,76.03,11.24,76.91},{23.79,83.81,24.64,84.59},{24.76,72.78,26.46,74.4},{27.85,77.07,28.27,77.55},{22.28,73.35,23.46,74.02},{30.47,76.78,30.93,77.18},{29.16,76.64,29.49,77.17},{23.82,79.74,25.13,80.68},{26.94,93.21,27.67,94.22},{18.75,76.21,19.87,77.12},{26.58,83.83,27.52,84.76},{21.76,86.56,22.95,87.89},{21.97,84.99,22.88,86.05},{23.4,71.03,24.14,72.48},{9.07,76.48,9.49,77.28},{29.79,75.93,30.68,76.8},{25.2,84.69,25.73,86.07},{11.05,78.63,11.52,79.17},{25.2,93.33,25.7,93.95},{25.45,94.19,25.92,94.9},{28.12,79.6,28.89,80.44},{29.44,79.82,30.81,81.04},{21.23,69.37,21.98,70.15},{14.95,78.74,16.32,80.48},{23.53,74.3,24.51,74.99},{25.57,81.33,26.18,82.44},{9.85,78.45,10.73,79.26},{33.78,74.78,34.08,75.19},{33.38,73.88,34.01,74.56},{17.89,73.32,19.39,75.16},{26.26,84.49,27.02,85.3},{21.6,87.43,22.51,88.19},{22.21,86.07,23.02,86.89},{19.46,85.14,20.2,86.38},{25.43,86.99,26.13,87.87},{22.7,85.83,23.7,86.9},{25.81,80.68,26.6,81.62},{15.55,76.24,16.56,77.6},{21.34,82.93,22.79,83.81},{17.85,72.84,19.13,73.67},{19.78,81.53,21.88,82.98},{22.79,77.36,23.75,78.82},{33.06,74.13,33.57,74.66},{23.46,76.2,24.28,77.23},{21.53,70.03,23.18,71.51},{20.11,80.39,21.83,81.22},{24.74,73.5,26.02,74.39},{12.24,77.07,13.19,77.64},{9.09,78.22,9.95,79.19},{33.1,75.01,33.58,75.42},{23.42,85.2,23.94,85.89},{28.38,78.84,29.17,79.41},{22.88,84.87,23.72,85.9},{16.85,77.36,17.71,78.86},{23.09,74.52,23.92,75.69},{16.49,73.03,18.07,73.87},{18.91,82.88,19.97,84.03},{32.87,74.52,33.53,75.08},{24.32,81.05,25.19,82.31},{27.97,76.28,28.47,76.86},{25.64,91.34,26.12,92.3},{28.68,76.21,29.1,76.9},{24.5,83.5,25.37,84.47},{30.22,78.82,30.81,79.36},{30.75,76.29,31.44,76.74},{23.06,72.73,24.49,73.66},{23.17,78.06,24.45,79.35},{29.56,77.12,30.41,77.94},{25.6,86.31,26.08,86.88},{24.72,87.45,25.35,87.98},{30.35,76.52,30.94,76.94},{21.94,92.82,22.81,93.21},{11.31,77.65,11.97,78.84},{25.46,85.53,26.09,86.42},{32.4,74.83,32.77,75.37},{20.92,83.8,22.19,84.77},{16.71,73.69,17.63,75.68},{29.73,75.56,30.69,76.2},{26.41,82.84,27.11,83.23},{25.18,82.19,25.54,82.71},{22.48,85.51,23.15,86.25},{25.62,84.41,26.22,85.2},{17.09,73.53,18.18,74.91},{23.96,80.36,25.19,81.39},{22.56,76.45,23.69,78.03},{24.59,93.67,25.62,94.49},{21.59,79.2,22.95,80.29},{23.0,92.68,23.52,93.2},{23.02,81.0,24.28,81.98},{30.97,75.78,31.28,76.52},{27.47,79.33,28.46,80.37},{23.11,75.69,24.33,77.05},{24.97,85.6,25.28,85.99},{26.4,85.18,26.65,85.39},{25.27,76.48,26.22,77.67},{30.76,76.99,31.71,78.31},{13.46,74.63,14.65,75.88},{24.84,77.01,25.92,78.47},{27.27,81.62,27.95,82.21},{33.62,74.56,33.83,75.03},{26.71,94.4,27.27,95.37},{27.02,82.44,27.5,83.29},{23.81,81.31,24.62,82.38},{27.12,74.68,28.2,76.1},{22.34,84.01,22.84,85.08},{15.61,73.31,16.66,74.21},{23.78,81.88,24.7,82.82},{30.38,77.02,31.02,77.83},{24.33,72.25,25.29,73.16},{29.22,74.47,29.99,75.3},{26.27,85.24,26.87,85.83},{27.1,80.3,27.91,81.42},{9.52,78.13,10.41,79.01},{25.89,84.01,26.37,84.78},{30.75,76.6,31.37,77.25},{17.11,74.61,18.55,76.43},{20.53,83.45,21.18,84.28},{23.87,82.53,24.91,83.55},{28.81,76.48,29.29,77.23},{26.5,92.33,27.04,93.79},{28.4,77.11,28.61,77.35},{27.08,88.27,27.53,88.54},{21.54,88.01,22.63,88.82},{25.14,90.23,25.55,90.97},{14.75,73.78,15.49,74.34},{22.94,91.32,23.76,91.87},{28.48,76.84,28.67,77.21},{18.08,83.41,19.17,84.77},{33.98,74.67,34.35,75.19},{25.97,81.54,26.66,82.69},{21.59,83.54,22.53,85.38},{26.0,86.4,26.56,87.11},{20.81,72.58,21.56,73.7},{22.12,70.95,23.53,72.19},{22.63,82.5,24.11,84.08},{24.49,93.15,25.46,93.96},{31.09,74.51,31.57,75.28},{27.46,91.55,27.86,92.48},{30.05,77.94,30.88,79.03},{18.99,72.64,20.23,73.8},{10.14,78.78,11.18,79.57},{20.56,73.47,21.08,73.94},{11.19,76.23,11.71,77.02},{9.52,77.17,10.22,77.73},{12.95,79.29,13.56,80.35},{8.29,76.67,8.86,77.28},{10.27,79.27,11.02,79.78},{8.32,77.67,9.37,78.39},{24.24,93.84,24.73,94.15},{10.17,75.95,10.78,76.9},{24.44,78.43,25.56,79.35},{27.22,95.22,27.97,96.02},{26.65,95.18,27.26,95.7},{10.29,78.16,11.4,79.01},{8.14,77.15,9.42,77.98},{10.23,77.05,11.35,77.89},{11.98,78.63,12.87,79.76},{25.68,75.11,26.56,76.33},{25.87,94.56,26.48,95.18},{12.75,76.35,14.34,77.52},{23.81,73.01,25.12,74.44},{26.48,91.71,26.93,92.41},{32.63,74.95,33.16,75.82},{28.71,78.73,29.38,80.16},{12.98,74.58,13.98,75.2},{22.83,75.13,23.76,76.25},{24.48,94.12,25.69,94.75},{23.19,80.53,24.09,81.35},{31.3,75.94,31.86,76.48},{26.11,80.05,27.03,81.05},{28.15,94.19,29.35,95.41},{27.76,93.19,28.71,94.6},{25.24,87.81,26.5,88.53},{13.93,74.05,15.53,75.11},{30.47,77.75,31.47,79.41},{21.82,72.85,22.81,74.29},{25.48,85.07,26.02,85.64},{20.12,72.73,20.74,73.49},{25.17,82.67,25.58,83.19},{12.26,78.41,13.21,79.79},{23.35,77.26,24.35,78.3},{11.5,78.63,12.45,80.0},{9.08,77.34,9.78,79.49},{17.25,81.87,18.55,83.49},{17.83,83.0,19.16,83.82},{17.32,78.83,18.61,80.7},{20.29,78.06,21.36,79.22},{19.85,76.61,20.76,77.69},{11.45,75.78,11.98,76.44},{28.61,76.95,28.7,77.19},{27.11,88.02,27.62,88.36},{25.2,89.82,25.95,90.44},{16.3,80.86,17.5,81.87},{26.89,92.02,27.78,92.91},{25.17,90.75,25.86,91.83},{21.37,75.21,22.55,76.24},{27.55,93.97,29.04,94.96},{23.27,91.15,24.23,91.79},{25.92,93.95,26.56,94.39},{16.19,76.29,16.95,77.48},{19.43,77.29,20.7,79.15},{25.76,94.34,26.29,94.72},{21.7,80.83,22.52,81.56},{29.93,77.07,30.5,77.61},{29.06,75.94,29.84,76.76},{26.38,76.15,27.23,77.08},{26.03,76.49,27.0,77.39},{25.74,75.98,26.72,76.98},{25.05,84.45,25.32,84.88},{24.99,84.86,25.32,85.22},{10.27,79.51,11.42,79.88},{10.83,79.72,11.0,79.85},{11.77,79.6,12.05,79.88},{16.7,82.18,16.76,82.31},{20.81,73.2,21.57,74.33},{6.76,92.72,9.26,93.95},{10.51,92.21,12.32,93.88},{12.06,92.65,13.68,94.28},{8.26,71.73,12.4,74.12} };

int main()
{
    freopen("results_kolkata.csv","w",stdout);
    freopen("test_gurgaon.txt","w",stderr);
    // int no_of_file = 3;

    // minlat = latlon[106][0];
    // maxlat = latlon[106][2];
    // minlon = latlon[106][1];
    // maxlon = latlon[106][3];

    minlat = 22.4900;
    maxlat = 22.6300;
    minlon = 88.2700;
    maxlon = 88.4100;

    cols = roundl((long double)(maxlon-minlon)*100.0);
    rows = roundl((long double)(maxlat-minlat)*100.0);
    // cerr<<"cols:"<<cols<<endl;
    // cerr<<"rows:"<<rows<<endl;
    long long tot_grids = roundl(((long double)(maxlat-minlat)*100.0)*((long double)(maxlon-minlon)*100.0));
    // cout<<"total grids:"<<tot_grids<<endl;
    // cout<<"here"<<endl;
    string inp_file ="",out_file ="";
    inp_file+="kolkata_processed";
    // inp_file+=to_string(i);
    inp_file+=".osm";
    // out_file+="dataout";
    // out_file+=to_string(i);
    // out_file+=".txt";
    ifstream inp_handle;
    inp_handle.open(inp_file.c_str());
    //freopen(inp_file.c_str(),"r",stdin);

    //freopen("tess.txt","w",stdout);
    // cout<<i<<" ";
    auto inp_time_start = high_resolution_clock::now();
    cout<<fixed;
    cout<<setprecision(4);
    //cout<<inp_file<<endl;
    while(getline(inp_handle,s))
    {
        s=ignore_leading_space(s);
        // cout<<s<<endl;
        if(s[s.size()-1]=='>')
        {
            s = s.substr(0,s.size()-1);
            if(s[s.size()-1]=='/')
            {
                s = s.substr(0,s.size()-1);

            }
            s+=' ';
            s+='>';
        }
        if(s.size()==0) continue;
        if(s.substr(1,4) =="node")
            node();
        else if(s.substr(1,3) =="way")
            way(inp_handle);
        else
            continue;

    }
    cerr<<"Nodes:\n";
    for(auto u:node_store){
        cerr<<u["id"]<<" ";
    }
    cerr<<'\n';
    for(int i=0;i<tot_grids;i++){

        cerr<<"No "<<i<<'\n';
        for(auto u:node_data[i])
            cerr<<u["id"]<<" ";
        cerr<<'\n';
        for(auto e:edge_list_all[i]){
            cerr<<e.first.first+1<<" "<<e.first.second+1<<'\n';
        }
        cerr<<'\n';
    }
    cerr<<"Total grids:"<<tot_grids<<'\n';
    cerr<<"Cols:"<<cols<<'\n';
    cerr<<"Rows:"<<rows<<'\n';

    cout<<"S.No,ThreeWay,FourWay,MoreThanFourWay,RoadLength,WalkabilityRatio\n";
    // cout<<"input complete!"<<endl;
    // cerr<<"node store size:"<<node_store.size()<<endl;
    for(long long i=0;i<tot_grids;i++){
        cout<<i<<",";
        long long iter=0;
        set<pair<long long, long long>> temp_edge_set;
        add_edges(i,temp_edge_set,iter);
        if(i==10){
            cerr<<"set:\n";
            for(auto u:temp_edge_set){
                cerr<<u.first<<" "<<u.second<<'\n';
            }
        }
        find_intersections();
        total_road_length(i);
        if((i+1)<tot_grids)
            add_edges(i+1,temp_edge_set,iter);
        if(i==10){
            cerr<<"set:\n";
            for(auto u:temp_edge_set){
                cerr<<u.first<<" "<<u.second<<'\n';
            }
        }
        if((i-1)>=0)
            add_edges(i-1,temp_edge_set,iter);
        if(i==10){
            cerr<<"set:\n";
            for(auto u:temp_edge_set){
                cerr<<u.first<<" "<<u.second<<'\n';
            }
        }
        if((i-rows)>=0)
            add_edges(i-rows,temp_edge_set,iter);
        if(i==10){
            cerr<<"set:\n";
            for(auto u:temp_edge_set){
                cerr<<u.first<<" "<<u.second<<'\n';
            }
        }
        if(((i-rows+1)>=0) and ((i-rows+1)<tot_grids))
            add_edges(i-rows+1,temp_edge_set,iter);
        if(i==10){
            cerr<<"set:\n";
            for(auto u:temp_edge_set){
                cerr<<u.first<<" "<<u.second<<'\n';
            }
        }
        if(((i-rows-1)>=0) and ((i-rows-1)<tot_grids))
            add_edges(i-rows-1,temp_edge_set,iter);
        if(i==10){
            cerr<<"set:\n";
            for(auto u:temp_edge_set){
                cerr<<u.first<<" "<<u.second<<'\n';
            }
        }
        if((i+rows)<tot_grids)
            add_edges(i+rows,temp_edge_set,iter);
        if(i==10){
            cerr<<"set:\n";
            for(auto u:temp_edge_set){
                cerr<<u.first<<" "<<u.second<<'\n';
            }
        }
        if((i+rows+1)<tot_grids)
            add_edges(i+rows+1,temp_edge_set,iter);
        if(i==10){
            cerr<<"set:\n";
            for(auto u:temp_edge_set){
                cerr<<u.first<<" "<<u.second<<'\n';
            }
        }
        if(((i+rows-1)<tot_grids) and ((i+rows-1)>=0))
            add_edges(i+rows-1,temp_edge_set,iter);
        if(i==10){
            cerr<<"set:\n";
            for(auto u:temp_edge_set){
                cerr<<u.first<<" "<<u.second<<'\n';
            }
        }
        if(i==10){
            cerr<<"check:\n";
            for(auto u:current_map){
                cerr<<node_store[u.first]["id"]<<'\n';
            }
        }
        walkability_ratio_v2(i);
        
        refresh_graph(i);
        if((i+1)<tot_grids)
            refresh_graph(i+1);
        if((i-1)>=0)
            refresh_graph(i-1);

        if((i-rows)>=0)
            refresh_graph(i-rows);
        if(((i-rows+1)>=0) and ((i-rows+1)<tot_grids))
            refresh_graph(i-rows+1);
        if(((i-rows-1)>=0) and ((i-rows-1)<tot_grids))
            refresh_graph(i-rows-1);
        
        if((i+rows)<tot_grids)
            refresh_graph(i+rows);
        if((i+rows+1)<tot_grids)
            refresh_graph(i+rows+1);
        if(((i+rows-1)<tot_grids) and ((i+rows-1)>=0))
            refresh_graph(i+rows-1);
        current_map.clear();
    }
    long double overall_road_length = 0.0;
    for(auto u:set_of_all_edges){
        overall_road_length+=(find_distance(find_lat(u.first),find_lon(u.first),find_lat(u.second),find_lon(u.second)))/2.0;
    }
    cerr<<fixed<<setprecision(10);
    cerr<<"overall_road_length:"<<overall_road_length<<'\n';
    for(int i=0;i<5000;i++){
        if(revisions[i]){
            cerr<<i<< " "<<revisions[i]<<'\n';
        }
    }
        //find_blocks();
        // time printing
        /*
        auto inp_time_end = high_resolution_clock::now();
        cout<<"\nInput Time:"<<duration_cast<microseconds>(inp_time_end-inp_time_start).count()<<" microseconds\n";

        auto rabt_time_start = high_resolution_clock::now();
        process_roundabouts();
        add_edges();
        auto rabt_time_end = high_resolution_clock::now();
        cout<<"\nGraph Creation Time:"<<duration_cast<microseconds>(rabt_time_end-rabt_time_start).count()<<" microseconds\n";
        auto int_time_start = high_resolution_clock::now();
        find_intersections();
        auto int_time_end = high_resolution_clock::now();
        cout<<"\nIntersections Finding Time:"<<duration_cast<microseconds>(int_time_end-int_time_start).count()<<" microseconds\n";
        //identify_blocks();
        auto walk_time_start = high_resolution_clock::now();
        walkability_ratio();
        auto walk_time_end = high_resolution_clock::now();
        cout<<"\nWalkability Ratio Time:"<<duration_cast<microseconds>(walk_time_end-walk_time_start).count()<<" microseconds\n";

        auto road_time_start = high_resolution_clock::now();
        total_road_length();
        auto road_time_end = high_resolution_clock::now();
        cout<<"\nRoad Length Time:"<<duration_cast<microseconds>(road_time_end-road_time_start).count()<<" microseconds\n";

        auto block_time_start = high_resolution_clock::now();
        find_blocks();
        auto block_time_end = high_resolution_clock::now();
        cout<<"\nBlocks Finding Time:"<<duration_cast<microseconds>(block_time_end-block_time_start).count()<<" microseconds\n";

        auto end = high_resolution_clock::now();
        auto time_taken = duration_cast<microseconds>(end-inp_time_start);
        cout<<"\nTotal Time taken: "<<time_taken.count()<<" microseconds\n";
        */
}
