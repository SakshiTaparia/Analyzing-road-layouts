#include<bits/stdc++.h>
using namespace std;

/** A function to find the distance between two points given their latitude and longitude. We use the Haverline Formula **/
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


/**

NOTE:
UPDATE THESE LAT LONG VALUES AND FILE NAMES ONLY!


Districts

1. Kolkata, West Bengal: MaxLon:88.41, MinLon:88.27, MaxLat:22.63, MinLat:22.5 (314)
2. Chennai, Tamil Nadu: MaxLon:80.31, MinLon:80.2, MaxLat:13.14, MinLat:12.99 (111)
3. Bangalore, Karnataka: MaxLon:77.84, MinLon:77.33, MaxLat:13.23, MinLat:12.66 (47)
4. Mumbai, Maharashtra: MaxLon:72.98, MinLon:72.78, MaxLat:19.27, MinLat:18.85 (382 + 383)
5. Gurgaon, Haryana: MaxLon:77.24, MinLon:76.66, MaxLat:28.54, MinLat:28.2 (213)
6. Chandigarh, Chandigarh: MaxLon:76.85, MinLon:76.7, MaxLat:30.8, MinLat:30.67 (107)
7. Hyderabad, Andhra Pradesh: MaxLon:78.54, MinLon:78.39, MaxLat:17.48, MinLat:17.3 (232)
8. Pune, Maharashtra: MaxLon:75.16, MinLon:73.32, MaxLat:19.39, MinLat:17.89 (452)
9. Delhi, Delhi: MaxLon:77.35, MinLon:76.84, MaxLat:28.88, MinLat:28.4 (99 + 166 + 410 + 412 + 416 + 419 + 543 + 549 + 612)

**/

long double pixel_threshold = 0.5, grid_threshold = 0.5, walking_distance_circle=564.0;
long double lat1=12.648253440856934,lon1=77.32221984863281,lat2=13.223280906677246,lon2=77.82559967041016;

vector<vector<int>> find_urban_extent(vector<vector<int>> &bu_matrix)
{
	vector<vector<int>> label = bu_matrix;	//label to indicate the type of the current pixel
	int N = bu_matrix.size(), M = bu_matrix[0].size();
	// walking distance circle's radius is taken to be 564 metres
	// long double walking_distance_circle = 1000.0;
	int ri = floor((long double)((walking_distance_circle*N)/find_distance(lat1,lon1,lat2,lon1))), rj = floor((long double)((walking_distance_circle*M)/find_distance(lat1,lon1,lat1,lon2)));
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			if(bu_matrix[i][j]==0){
				label[i][j] = 3;
				continue;
			}
			int built = 0, tot = 0;
			for(int u=i-ri ; u<=i+ri;u++)
			{
				for(int v = j-rj; v<=j+rj; v++)
				{
					if(u<0 or u>=N or v<0 or v>=M)
						continue;
					if((rj*rj*((double)(u-i)*(u-i)) + ri*ri*((double)(v-j)*(v-j)) - ri*ri*rj*rj)<=0)
					{
						built+=(bu_matrix[u][v]==1);
					}
					tot++;
				}
			}
			//cout<<built<<" "<<tot<<'\n';
			long double perc = ((long double)(built))/tot;
			if(perc>=pixel_threshold)
				label[i][j] = 0;	//urban
			else if(perc>=(pixel_threshold/2.0))
				label[i][j] = 1;	//suburban
			else
				label[i][j] = 2;	//rural
		}
	}
	//labelled all pixels
	return label;
}

void preprocess(vector<vector<int>> &matrix,bool year2019){
	if(year2019){
		for(auto &u:matrix){
			for(auto &v:u){
				if(v==255)
					v=2;
				else if((v==200) or (v==100))
					v=1;
				else
					v=0;
			}
		}
	}
	else{
		for(auto &u:matrix){
			for(auto &v:u){
				if((v==255) or (v==100))
					v=2;
				else if(v==200)
					v=1;
				else
					v=0;
			}
		}
	}
}

vector<vector<int>> rotate_matrix(vector<vector<int>> &matrix){
	int new_i=0,new_j=0;
	int n = matrix.size(),m=matrix[0].size();
	vector<vector<int>> new_mat(m,vector<int>(n));
	for(int j=0;j<m;j++){
		new_j=0;
		for(int i=n-1;i>=0;i--){
			new_mat[new_i][new_j++] = matrix[i][j];
			new_j%=n;
		}
		new_i++;
	}
	return new_mat;
}



int main()
{
	freopen("new_input.txt","r",stdin);
	freopen("bangalore_testing_2014.txt","w",stdout);
	vector<vector<int>> bu_matrix;
	string cur_line;
	deque<deque<int>> inp_matrix;
	while(getline(cin,cur_line)){
		deque<int> temp;
		int i=0,j=0;
		while(i<cur_line.size()){
			int val = 0;
			j = i;
			while((j<cur_line.size()) and (cur_line[j]!=' ')){
				val*=10;
				val+=(cur_line[j]-'0');
				j++;
			}
			temp.push_back(val);
			i = j+1;
		}	
		inp_matrix.push_back(temp);
	}
	long double rounded_lat1,rounded_lon1,rounded_lat2,rounded_lon2;
	rounded_lat1 = floor(100.0*lat1)/100.0;
	rounded_lon1 = floor(100.0*lon1)/100.0;
	rounded_lat2 = ceil(100.0*lat2)/100.0;
	rounded_lon2 = ceil(100.0*lon2)/100.0;

	int top,bottom,left,right;

	top = ((long double)inp_matrix.size())*((rounded_lat2-lat2)/(lat2-lat1));
	bottom = ((long double)inp_matrix.size())*((lat1-rounded_lat1)/(lat2-lat1));
	left = ((long double)inp_matrix[0].size())*((lon1-rounded_lon1)/(lon2-lon1));
	right = ((long double)inp_matrix[0].size())*((rounded_lon2-lon2)/(lon2-lon1));

	int initial_length = inp_matrix[0].size();

	for(int i=0;i<top;i++){
		deque<int> q;
		for(int j=0;j<initial_length;j++)
			q.push_back(0);
		inp_matrix.push_front(q);
	}
	for(int i=0;i<bottom;i++){
		deque<int> q;
		for(int j=0;j<initial_length;j++)
			q.push_back(0);
		inp_matrix.push_back(q);
	}
	for(int i=0;i<inp_matrix.size();i++){
		for(int j=0;j<left;j++)
			inp_matrix[i].push_front(0);
		for(int j=0;j<right;j++)
			inp_matrix[i].push_back(0);
	}
	for(int i=0;i<inp_matrix.size();i++){
		vector<int> temp;
		for(int j=0;j<inp_matrix[0].size();j++){
			temp.push_back(inp_matrix[i][j]);
		}
		bu_matrix.push_back(temp);
	}
	int N = bu_matrix.size(), M = bu_matrix[0].size();
	lat1 = rounded_lat1;
	lon1 = rounded_lon1;
	lat2 = rounded_lat2;
	lon2 = rounded_lon2;

	for(auto u:bu_matrix)
	{
		for(int i=0;i<u.size();i++)
			{
				cout<<u[i];
				if(i<u.size()-1)cout<<" ";
			}
		cout<<'\n';
	}
	preprocess(bu_matrix,false);
	

	vector<vector<int>> out_mat = find_urban_extent(bu_matrix);
	// vector<vector<int>> out_mat = (bu_matrix);
	freopen("bangalore_urban_extent_2014.txt","w",stdout);
	for(auto u:out_mat)
	{
		for(int i=0;i<u.size();i++)
			{
				cout<<u[i];
				if(i<u.size()-1)cout<<" ";
			}
		cout<<'\n';
	}
	
	long long height, width;
	height = llroundl(abs(lat2-lat1)*100.0);
	width = llroundl(abs(lon2-lon1)*100.0);
	swap(height,width);
	swap(N,M);
	out_mat = rotate_matrix(out_mat);
	// cerr<<"rot:\n";
	// for(auto u:out_mat){
	// 	for(auto v:u){
	// 		cerr<<v<<" ";
	// 	}
	// 	cerr<<'\n';
	// }
	cout<<fixed;
	cout<<setprecision(5);
	set<int> use,sub,urb,rur,discard;
	set<int> set_of_grids, discarded_grids;
	map<int,long double> urb_perc,sub_perc;
	for(long double i=0;i<height;i+=1.0)
	{
		for(long double j=0;j<width;j+=1.0)
		{
			int tot = 0, cnt=0, cnt1=0, cnt2=0, cnt3=0;
			//cout<<floor((1858.0*i)/50.0)<<" "<<floor((1858.0*(i+1.0))/50.0)<<" "<<floor((2228.0*j)/60.0)<<" "<<floor((2228.0*(j+1.0))/60.0)<<'\n';
			for(int u=floor((N*i)/height) ; u<floor((N*(i+1.0))/height);u++)
			{
				for(int v=floor((M*j)/width);v<floor((M*(j+1.0))/width);v++)
				{
					if((u>=0 and u<out_mat.size()) and (v>=0 and v<out_mat[0].size()))
					{
						cnt+=((out_mat[u][v]==0) or (out_mat[u][v]==1));
						cnt1+=(out_mat[u][v]==0);
						cnt2+=(out_mat[u][v]==1);
						cnt3+=(out_mat[u][v]==3);
						tot++;
					}
				}
			}
			long double perc = ((long double)cnt)/tot;
			long double include_perc = ((long double)cnt3)/tot;
			int grid_no = roundl((i)*(width) + (j) );
			//cout<<(i-0.5)*59.0 + (j-0.5) + 1.0<<" "<<perc<<endl;
			if(perc>=grid_threshold)
			{
				// use.push_back((int)((i-0.5)*(width-1) + (j-0.5) + 1));
				use.insert(grid_no);
				if(cnt1>cnt2){
					urb.insert(grid_no);
				}
				else{
					sub.insert(roundl((i)*(width) + (j) ));
				}
				urb_perc[grid_no] = ((long double)cnt1)/tot;
				sub_perc[grid_no] = ((long double)cnt2)/tot;
			}
			else if(include_perc<0.5){
				urb_perc[grid_no] = ((long double)cnt1)/tot;
				sub_perc[grid_no] = ((long double)cnt2)/tot;
				rur.insert(roundl((i)*(width) + (j) ));
			}
			else{
				discard.insert(roundl((i)*(width) + (j) ));
				discarded_grids.insert(roundl((i)*(width) + (j) ));
			}
			set_of_grids.insert(grid_no);

		}
	}

	freopen("usable_grids_bangalore_2014.txt","w",stderr);

	// cerr<<set_of_grids.size()<<endl;
	// for(auto u:set_of_grids){
	// 	cerr<<u<<'\n';
	// }
	cerr<<use.size()<<'\n';
	for(auto u:use)
	{
		cerr<<u<<'\n';
	}
	cerr<<urb.size()<<'\n';
	for(auto u:urb)
	{
		cerr<<u<<'\n';
	}
	cerr<<sub.size()<<'\n';
	for(auto u:sub)
	{
		cerr<<u<<'\n';
	}
	cerr<<rur.size()<<'\n';
	for(auto u:rur){
		cerr<<u<<'\n';
	}
	cerr<<discard.size()<<'\n';
	for(auto u:discard){
		cerr<<u<<'\n';
	}

	cerr<<"suburban/urban:"<<((long double)sub.size())/urb.size()<<'\n';
	cerr<<"rural/(urban+suburban):"<<((long double)rur.size())/(sub.size()+urb.size())<<'\n';
	freopen("urb_sub_perc_bangalore_2014.csv","w",stdout);
	cout<<"grid_no,urb_perc,sub_perc,district,s/r\n";
	for(auto u:set_of_grids){
		cout<<u<<',';
		if(discarded_grids.find(u)==discarded_grids.end()){
			cout<<urb_perc[u]<<","<<sub_perc[u]<<',';
		}
		else{
			cout<<"NA,NA,";
		}
		cout<<"bangalore_2014,";
		if(use.find(u)!=use.end()){
			cout<<"Selected\n";
		}
		else
			cout<<"Rejected\n";
	}
}