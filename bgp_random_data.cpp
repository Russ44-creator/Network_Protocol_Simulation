#include<iostream>
#include<fstream>
#include<algorithm>
#include<string>
#include<vector>
#include<cstdio>
#include<unistd.h>
#include<sys/time.h>
#include <ctime> 
#include <chrono>
#include <cstring>
#include <queue>
#include <climits>


using namespace std;

#define INF INT_MAX //Infinity


// routing table
struct Routing{
	string name; //Router Name
	int N;  // Destination network N
	int d; // Distance
	string X;  //Next Router
    int as; 
};

const int sz=10001; //Maximum possible number of vertices. Preallocating space for DataStructures accordingly
vector<Routing> routers[1000000];  // route table
char idToRname[] = "ABCDEF";
//vector< vector<int> > ve(10000, vector<int>(10000, INF));
vector<pair<int,int> > ve[sz]; 
//int ve[100000][100000] = {};  // network graph

string num_to_name(int source, int k){
		string name_router;
		int layer_2 = k * k / 2;
		if(source < layer_2){
			auto temp = std::to_string(source);
			name_router = 'A' + temp;
		}
		else if(source >= layer_2 && source < layer_2 * 2){
			auto temp = std::to_string(source - layer_2);
			name_router = 'B' + temp;
		}
		else{
			auto temp = std::to_string(source - layer_2 * 2);
			name_router = 'C' + temp;
		}
		return name_router;
}

vector<Routing> Send2(vector<Routing> r, string c){
	for(int i = 0;i < r.size();i++){
		r[i].d++;
		r[i].X = c;
	}
	return r;
}


// init fattree
void fattree(int k){
	int core = (k / 2) ^ 2;
	int pod = k;
	int switch_in_pod = k;
	int server_each_pod = k * k /4;
	int layer_2 = k * k / 2;
	Routing rp;
	int N = core + pod * pod;
	for (int i = 0; i < layer_2; i++){
		auto temp = std::to_string(i);
		rp.name = 'A' + temp;
		for (int j = 0; j < k / 2; j++){
			rp.N = k / 2 * i + j;
			rp.d = 1;
			rp.X = '-';
            rp.as = -1;
			routers[i].push_back(rp);
		}	
	}
}
	
void fattree_graph(int k){
	int core = (k / 2) ^ 2;
	int pod = k;
	int layer_2 = k * k / 2;
	int switch_in_pod = k;
	int server_each_pod = k * k /4;

	// access layer
    //edgeport
	for (int i = 0; i < layer_2; i++){
		for (int j = 0; j < pod/2; j++){
			int distance = rand() % 10 + 1;
			//ve[i][layer_2 + j + (i / (k/2)) * (k/2)] = rand() % 10 + 1;
			//ve[layer_2 + j + (i / (k/2)) * (k/2)][i] = ve[i][layer_2 + j + (i / (k/2)) * (k/2)];
			ve[i].push_back(make_pair(layer_2 + j + (i / (k/2)) * (k/2), distance));
			ve[layer_2 + j + (i / (k/2)) * (k/2)].push_back(make_pair(i, distance));
		}		
	}
	// Convergence layer
	for (int i = layer_2; i < layer_2 * 2; i++){
		/*
		// access layer
		for(int j = 0; j < pod/2; j++){
			ve[i][(j + (i - layer_2)/(k/2) * (k/2))] = rand() % 10 + 1;
		}
		*/
		// core layer
		//for(int j = layer_2 * 2; j < layer_2*2 + core; j++){
		int group = (i - layer_2) % (k/2);
		for(int j = 0; j < k/2; j++){
			int distance = rand() % 10 + 1;
			//ve[i][(2*layer_2 + (k/2) * group + j)] = rand() % 10 + 1;
			//ve[(2*layer_2 + (k/2) * group + j)][i] = ve[i][(2*layer_2 + (k/2) * group + j)];
			ve[i].push_back(make_pair(2*layer_2 + (k/2) * group + j, distance));
			ve[2*layer_2 + (k/2) * group + j].push_back(make_pair(i, distance));
		}		
	}
	
}

string name_to_num(int num, int k){
	int layer_2 = k * k / 2;
	if(num < layer_2){
			auto temp = std::to_string(num);
			return 'A' + temp;
		}
	else if(num >= layer_2 && num < layer_2 * 2){
		auto temp = std::to_string(num - layer_2);
		return 'B' + temp;
	}
	else{
		auto temp = std::to_string(num - layer_2 * 2);
		return 'C' + temp;
	}

}

void printR2(int k){
	vector<Routing> r;
	for(int i = 0;i < k * k * 5 / 4;i++){
		r = routers[i];
		string temp = name_to_num(i, k);
		cout<<"***********************routing table "<<temp<<"************************\n";
		cout<<"To target N\t\tDistanced\t\tNext RouterX\t\tNext As\n";
		for(int i = 0;i < r.size();i++){
			cout<<r[i].N<<"\t\t\t"<<r[i].d<<"\t\t"<<r[i].as<<endl;
		}
	}
}

void printG(int k){
	cout<<"**********************Network Graph********************\n";
	for(int i = 0;i < k * k * 5 / 4;i++){
		cout << i << " ";
		for(int j = 0;j < ve[i].size();j++){
			if (ve[i][j].second != INF){
				cout << "router " << ve[i][j].first << " distance: " << ve[i][j].second << ", ";
				if(j < ve[i].size() - 1) cout<<" ";
			}
		}
		cout << endl;
	}
}


int dis[sz]; //Stores shortest distance
int next_node[sz][sz];
bool vis[sz]={0}; //Determines whether the node has been visited or not


void bgp(int k){

    int layer_2 = k * k / 2;
	int switch_in_pod = k;

    //edge port
    
	// select next node which is nearest
	// for the least distance
    int next_node = 0;
    for(int i = 0; i < layer_2; i++){
        int min = INF;
        for(int j = 0;j < ve[i].size();j++){
            if (ve[i][j].second < min){
                min = ve[i][j].second;
                next_node = ve[i][j].first;
            }
        }
       
        Routing rp;
		rp.name = num_to_name(i, k);
		rp.N = -99; //others
        rp.d = min;
        rp.X = next_node;
        rp.as = 10000 + layer_2 + i % (k/2);
        routers[i].push_back(rp);
    }

    Routing rp;

    // aggregation layer
    int min = INF;
    for(int i = layer_2; i <  2 * layer_2; i++){
            for(int j = 0;j < ve[i].size();j++){
                min = INF;
				
				//direct to the edgeports
                if (ve[i][j].first < layer_2){
                     int dis = ve[i][j].second;
                     int next_node = ve[i][j].first;
                     for(int t = 0; t < k / 2; t++){
                        rp.name = num_to_name(i, k);
                        rp.N = t + next_node * k / 2; 
                        rp.d = dis;
                        rp.X = num_to_name(next_node, k);
                        rp.as = 10000 + next_node;
                        routers[i].push_back(rp);
                     }

                }
				// to the core ports
				// choose the least distance

                if (ve[i][j].first > layer_2){
                    if (ve[i][j].second < min){
                        min = ve[i][j].second;
                        next_node = ve[i][j].first;
                    }
                }
           
        }
        
        
        rp.name = num_to_name(i, k);
        rp.N = -99; //others
        rp.d = min;
        rp.X = num_to_name(next_node, k);
        rp.as = 10000 + layer_2 + k;
        routers[i].push_back(rp);
    
    
    }

    // core layer
    for (int i = 2 * layer_2; i < 2 * layer_2 + k * k / 4; i++){
        for(int j = 0;j < ve[i].size();j++){

			//least distance
			// if same distance, choose the last find
            int dis = ve[i][j].second;
            next_node = ve[i][j].first;
            for (int t = 0; t < k * k / 4; t++){
                rp.name = num_to_name(i, k);
                rp.N = ((next_node - layer_2) / (k/2)) * (k / 2) + t; 
                rp.d = dis;
                rp.X = num_to_name(next_node, k);
                rp.as = 10000 + (next_node - layer_2) / (k/2);
                routers[i].push_back(rp);
            }
        }
    }
    
}

void Dijkstra(int source, int k){
	bool vis[sz]={0};
	int server_num = k * k / 2;
	for(int i = 0;i < sz;i++) //Set initial distances to Infinity
        dis[i] = INF;
	 //Custom Comparator for Determining priority for priority queue (shortest edge comes first)
    class prioritize{public: bool operator ()(pair<int, int>&p1 ,pair<int, int>&p2){return p1.second>p2.second;}};
    priority_queue<pair<int,int> ,vector<pair<int,int> >, prioritize> pq; //Priority queue to store vertex,weight pairs
	pq.push(make_pair(source,dis[source]=0)); //Pushing the source with distance from itself as 0
	int switch_num = k * k * 5 / 4;
	next_node[source][source] = -2;
	while(!pq.empty()){
		pair<int, int> curr=pq.top(); //Current vertex. The shortest distance for this has been found
        pq.pop();
		int cv = curr.first,cw = curr.second; //'cw' the final shortest distance for this vertex
        if(vis[cv]) //If the vertex is already visited, no point in exploring adjacent vertices
            continue;
        vis[cv] = true; 
		//next_node[source][cv] = ;
        for(int i = 0;i < ve[cv].size();i++) //Iterating through all adjacent vertices
            if(!vis[ve[cv][i].first] && ve[cv][i].second + cw < dis[ve[cv][i].first]) {
			//If this node is not visited and the current parent node distance+distance from there to this node is shorter than the initial distace set to this node, update it
                pq.push(make_pair(ve[cv][i].first,(dis[ve[cv][i].first]=ve[cv][i].second+cw))); //Set the new distance and add to priority queue
				if (cv == source)
					next_node[source][ve[cv][i].first] = ve[cv][i].first;
				else
					next_node[source][ve[cv][i].first] = next_node[source][cv];
			}
	}
}



int main(){
	int k = 20;
	fattree(k);
	//printR2(k);
	fattree_graph(k);
	//printG(k);
	int T=0;
	int count = 0;
	int switch_num = k * k * 5 / 4;
	int layer_2 = k * k / 2;
	auto start = std::chrono::system_clock::now();

	bgp(k);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> time_spent = end-start;
	
	//printR2(k);
	cout << "\n time spent in OSPF: " << time_spent.count() << endl;;
	return 0;
}
