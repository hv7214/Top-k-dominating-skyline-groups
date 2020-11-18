#include<bits/stdc++.h>
using namespace std;

class Point {
public:
    string name;
    vector<int> dims;
    Point(vector<int> dims, string name) : name(name){
        this->dims = dims;
    }

    /*
    returns:
        1 -> if current point dominates other point
        0 -> if no one dominates
        -1 -> if other point dominates current point
    */
    int dominates(Point& other) {
        int g = 0, l = 0, e = 0;
        int d = other.dims.size();
        for(int i = 0; i < d; i++) {
            if(dims[i] < other.dims[i]) l++;
            else if(dims[i] > other.dims[i]) g++; 
            else e++;
        }
        if(g >= 1 && (g + e) == d) return 1;
        else if(l >= 1 && (l + e) == d) return -1;
        return 0;
    }

    bool operator<(const Point& other) const {
        return name.compare(other.name);
    }
};

// Dataset
vector<Point> D = {Point({10, 0, 0}, "Q1"), Point({1,7,1}, "Q2"), Point({2,3,3}, "Q3"), Point({7,0,0}, "Q4"), Point({1,5,0}, "Q5"), Point({0,0,4}, "Q7"), Point({0,4,1}, "Q6"), Point({1,1,2}, "Q8"), Point({0,1,3}, "Q9"), Point({1,0,3}, "Q10")};

class Group {
public:
    set<Point> points;
    int score;

    Group(set<Point> points) : points(points) {
        score = 0;
    }

    Group() {
        points = {};
    }

    bool operator==(const Group& other) const {
        if(other.points.size() != points.size()) return false;
        map<string, bool> hmap;
        
        for(auto p : points) {
            hmap[p.name] = true;
        }
        for(auto p : other.points) {
            if(hmap[p.name] == false) return false;
        }
        
        return true; 
    }
};

// Returns the score of the given group 
int getScore(Group G, vector<Point>& D) {
    map<string, bool> hmap;
    int score = 0;
    for(auto q : G.points) {
        for(auto p : D) {
            if(q.dominates(p) == 1 && !hmap[p.name]) {
                hmap[p.name] = true;
                score++;
            }
        }
    }
    return score;
}

// compartor for Group class
struct comp{ 
    bool operator()(const Group& curr, const Group& other) {
        return getScore(curr, D) > getScore(other, D);
    }
};

// Merge units, used in Unit+ algorithm
Group mergeUnits(Group g1, Group g2) {
    for(auto p : g2.points) {
        g1.points.insert(p);
    }
    return g1;
}

// prunes the input using Unit+ algorithm
vector<Point> inputPruning(vector<Point> D, int l) {
    int n = D.size();
    int m = n / l;
    vector<Point> D_pruned;
    vector<Group> units, groups;
    map<string, bool> hmap;

    // form units of every point in the dataset
    for(int i = 0; i < n; i++) {
        auto p = D[i];
        Group unit;
        for(int j = 0; j < n; j++) {
            auto q = D[j];
            if(i == j) unit.points.insert(p);
            else {
                if(q.dominates(p) == 1) {
                    unit.points.insert(q);
                }
            }
        }
        // if a point is dominated by more than l - 1 points,
        // then the point will not be in the dataset
        if(unit.points.size() <= l) D_pruned.push_back(p);
    }

    return D_pruned;
} 

// Returns the skyline(points which are dominated by any point in the dataset) of the dataset
vector<Point> findSkyline(vector<Point>& D) {
    int n = D.size();
    if(n == 0) return {};

    vector<Point> skyline = {D[0]};
    for(int i = 1; i < n; i++) {
        bool flag = true;
        for(int j = 0; j < skyline.size(); j++) {
            int res = skyline[j].dominates(D[i]);
            if(res == -1) {
                skyline.erase(begin(skyline) + j);
                j--;
            } 
            else if(res == 1) {
                flag = false;
            }
        }
        if(flag) {
            skyline.push_back(D[i]);
        }
    }
    return skyline;
}

// removes the skyline points from the dataset and returns residual
vector<Point> findResidual(vector<Point>& D, vector<Point>& skyline) {
    map<string, bool> hmap;
    for(auto s : skyline) hmap[s.name] = true;

    vector<Point> residual;
    for(auto p : D) {
        if(!hmap[p.name]) {
            residual.push_back(p);
        }
    }
    return residual;
}

// generate skyline candidate groups, used in TKD_permutation
void generateGroups(int pos, int l, int k, 
            priority_queue<Group, vector<Group>, comp>& PQ, 
            vector<Point> skyline, Group candidate) {
    if(l == 0) {
        return;
    }
    for(int i = pos; i < skyline.size(); i++) {
        auto temp = candidate;
        temp.points.insert(skyline[i]);

        if(PQ.size() < k) {
            PQ.push(temp);
        } else {
            if(getScore(temp, D) > getScore(PQ.top(), D)) {
                PQ.pop();
                PQ.push(temp);
            }
        }
        generateGroups(i+1, l-1, k, PQ, skyline, temp);
    }
}

void revive(int pos, int l, int k, priority_queue<Group, 
            vector<Group>, comp>& PQ, vector<Point> residual, Group candidate) {
    if(l == 0) {
        if(PQ.size() < k) {
            PQ.push(candidate);
        }
        else if(getScore(candidate, D) > getScore(PQ.top(), D)) {
            PQ.pop();
            PQ.push(candidate);
        }
    } else {
        for(int i = pos; i < residual.size(); i++) {
            auto temp = candidate;
            // check if temp points dominate residual[i]
            bool dom = true;
            for(auto p : temp.points) {
                int res = p.dominates(residual[i]);
                if(res == -1) {
                    dom = false;
                }
            }

            if(dom) {
                temp.points.insert(residual[i]);
                revive(i+1, l-1, k, PQ, residual, temp);
            }
        }
    }
}

vector<Group> TKD_permutation(vector<Point> D, int k, int l) {
    // prunes the dataset using unit algorithm
    D = inputPruning(D, l);
    // skyline points are those which are not dominated by any other point in dataset
    auto skyline = findSkyline(D);
    // subtract skyline from dataset(D - skyline)
    auto residual = findResidual(D, skyline);
    // helper data-structures
    vector<Group> candidate, bucket;
    // priority queue of size k which sorts groups on the basis of ascending scores 
    priority_queue<Group, vector<Group>, comp> PQ;
    // stores top-k dominating skyline groups
    vector<Group> topKSkylineGroups;
    // form candidate groups which can be in top-k dominating skyline groups
    generateGroups(0, l, k, PQ, skyline, Group(set<Point>()));

    while(PQ.size()) {
        auto G = PQ.top(); PQ.pop();

        if(G.points.size() < l) {
            candidate.push_back(G);
        } else {
            bucket.push_back(G);
        }
    }

    for(auto G : bucket) {
        PQ.push(G);
    }

    for(auto G : candidate) {
        if(PQ.size() < k || getScore(G, D) > getScore(PQ.top(), D)) {
            revive(0, l-G.points.size(), k, PQ, residual, G);
        }
    }

    while(PQ.size()) {
        topKSkylineGroups.push_back(PQ.top()); PQ.pop();
    }

    return topKSkylineGroups;
}

/* modify dimensions to improve prunning power of lemma 3 */
void modify_dims(vector<Point> &D){
    return;
}

bool comp_L1_norm(Point p1, Point p2){
    int lp1 = 0, lp2 = 0;
    for(auto d : p1.dims){
        lp1 += d;
    }
    for(auto d : p2.dims){
        lp2 += d;
    }
    return lp1 > lp2;
}

Point f_sum(Group grp){
    if(grp.points.size()==0){
        cout<<"f_sum error"<<endl;
        exit(1);
    }
    vector<int> new_dim;
    for(auto g : grp.points){
        if(new_dim.empty()){
            new_dim = g.dims;
        } else {
            for(int i=0;i<new_dim.size();i++){
                new_dim[i] += g.dims[i];
            }
        }
    }
    Point Q(new_dim, "Q");
    return Q;
}

vector<Group> skyline_groups(vector<Group> sky, vector<Group> grp){
    if(sky.size() == 0){
        return grp;
    }
    vector<Group> skyline;

    for(auto G : grp) {
        sky.push_back(G);
    }    
    for(auto G1 : sky) {
        bool flag = true;
        for(auto G2 : sky) {
            auto p1 = f_sum(G1);
            auto p2 = f_sum(G2);
            if(p2.dominates(p1) == 1) {
                flag = false;
                break;
            }   
        }
        if(flag) {
            skyline.push_back(G1);
        }
    }
    return skyline;
}

vector<Group> TKD_SUM(vector<Point> D, int k, int l) {
    D = inputPruning(D, l);
    modify_dims(D);
    sort(D.begin(), D.end(), comp_L1_norm);
    priority_queue<Group, vector<Group>, comp> PQ;
    vector<Group> Sky[D.size() + 1][l + 1];
    vector<Group> topKSkylineGroups;

    for(int j = 1;j<=D.size();j++){
        for(int i = min(j,l); i>=1;i--){
            if(i==1){
                Group g; g.points.insert(D[j-1]);
                vector<Group> dj; dj.push_back(g);
                Sky[j][1] = skyline_groups(Sky[j-1][1], dj);
            } else {
                vector<Group> temp;
                for(auto G : Sky[j-1][i-1]) {
                    G.points.insert(D[j-1]);
                    temp.push_back(G);
                }
                Sky[j][i] = skyline_groups(Sky[j-1][i], temp);
            }
        }
    }
    for(auto g : Sky[D.size()][l]){
        if(PQ.size() < k){
            PQ.push(g);
        } else {
            if(getScore(g, D) > getScore(PQ.top(), D)){
                PQ.pop();
                PQ.push(g);
            }
        }
    }
    while(PQ.size()) {
        topKSkylineGroups.push_back(PQ.top()); PQ.pop();
    }

    return topKSkylineGroups;
}

void findAllGroups(int pos, vector<Point> D, int groupSize, vector<Group>& groups, Group currGroup) {
    if(groupSize <= 0) {
        groups.push_back(currGroup);
        return;
    }

    for(int i = pos; i < D.size(); i++) {
        // include the current point
        auto include = currGroup; include.points.insert(D[i]);
        findAllGroups(pos+1, D, groupSize-1, groups, include);
        // exclude the current point
        findAllGroups(pos+1, D, groupSize, groups, currGroup);
    }
}

vector<int> findSkylineVector(const Group& G, int l) {
    vector<int> skylineVector;

    for(int i = 0; i < l; i++) {
        int p_max = INT_MIN;
        for(auto p : G.points) {
            p_max = max(p_max, p.dims[i]);
        }
        skylineVector.push_back(p_max);
    }
    return skylineVector;
}

bool sameSkylineVector(const Group& G1, const Group& G2, int d) {
    auto skylineVector1 = findSkylineVector(G1, d);
    auto skylineVector2 = findSkylineVector(G2, d);

    for(int i = 0; i < d; i++) {
        if(skylineVector1[i] != skylineVector2[i]) return false;
    }
    return true;
}

void getAllCombinations(vector<Group>& allCombinations, vector<vector<Point>> V, int pos, Group currGroup) {
    if(pos >= V.size()) {
        allCombinations.push_back(currGroup);
        return;
    }
    for(int i = 0; i < V[pos].size(); i++) {
        auto newGroup = currGroup; newGroup.points.insert(V[pos][i]);
        getAllCombinations(allCombinations, V, pos+1, newGroup);
    }
}

bool checkIfNotPresent(priority_queue<Group, vector<Group>, comp> PQ, Group G) {
    while(PQ.size()) {
        if(G == PQ.top()) return false;
        PQ.pop();
    }
    return true;
}

void constructGroups(Group G, int l, int k, priority_queue<Group, vector<Group>, comp>& PQ, vector<Point> D) {
    auto d = (*begin(G.points)).dims.size();
    auto v = findSkylineVector(G, d);   
    vector<vector<Point>> V(d);
    vector<Group> allCombinations;
    vector<Group> groups;

    for(int i = 0; i < d; i++) {
        for(auto p : D) {
            if(p.dims[i] == v[i]) {
                V[i].push_back(p);
            }
        }
    } 

    getAllCombinations(allCombinations, V, 0, Group());

    for(auto Gc : allCombinations) {
        if(Gc.points.size() <= l) {
            // find D - Gc
            auto temp = vector<Point>(begin(Gc.points), end(Gc.points));
            auto residual = findResidual(D, temp);
            // find all groups of size |l - Gc| from (D - Gc.points) set points
            findAllGroups(0, residual, l - Gc.points.size(), groups, Group());
            
            for(auto G2 : groups) {
                auto mergedGroup = mergeUnits(G2, Gc);

                if(PQ.size() < k && checkIfNotPresent(PQ, mergedGroup)) {
                    PQ.push(mergedGroup);
                }
                else if(getScore(mergedGroup, D) > getScore(PQ.top(), D) && checkIfNotPresent(PQ, mergedGroup)) {
                    PQ.pop();
                    PQ.push(mergedGroup);
                }
            }
        }
    }
}

vector<Group> TKD_MAX(vector<Point> D, int k, int l) {
	auto skyline = findSkyline(D);
    auto residual = findResidual(D, skyline);
    vector<Group> queue;
	vector<Group> topKSkylineGroups;
    vector<Group> groupCombinations;
	Group G(set<Point>(skyline.begin(), skyline.end()));
    vector<vector<vector<Group>>> Sky(skyline.size()+1, vector<vector<Group>>(l+1, vector<Group>()));

    // priority queue of size k which sorts groups on the basis of ascending scores 
    priority_queue<Group, vector<Group>, comp> PQ;

	if(G.points.size() <= l) {
        // find all groups of size l - |G| from residual points
        // NP-hard problem (Set-cover problem)
        findAllGroups(0, residual, l - G.points.size(), groupCombinations, Group());
        
        for(auto G1 : groupCombinations) {
            auto mergedGroup = mergeUnits(G, G1);
            topKSkylineGroups.push_back(mergedGroup);
            if(topKSkylineGroups.size() == k) {
                return topKSkylineGroups;
            }
        }
	}

    for(int j = 1; j <= skyline.size(); j++) {
        for(int i = 1; i <= min(j, l); i++) {
            if(i == 1) {
                Group g; g.points.insert(skyline[j-1]);
                vector<Group> tempGroup; tempGroup.push_back(g);
                Sky[j][i] = skyline_groups(Sky[j-1][i], tempGroup);
            } else {
                vector<Group> temp;
                for(auto G : Sky[j-1][i-1]) {
                    G.points.insert(skyline[j-1]);
                    temp.push_back(G);
                }
                Sky[j][i] = skyline_groups(Sky[j-1][i], temp);
            }
        }
    }
	
    // groups sorted in descending order
    auto groups_set = Sky[skyline.size()][l];

    // remove group from the groups set having same skyline/aggregate vector and less score
    for(auto Gi : groups_set) {
        bool takeFlag = true;
        for(auto Gj : groups_set) {
            if(Gi == Gj) continue;
            int score_Gi = getScore(Gi, D), score_Gj = getScore(Gj, D);
            if(score_Gi < score_Gj && sameSkylineVector(Gi, Gj, D[0].dims.size())) {
                takeFlag = false;
                break;
            }
        }
        if(takeFlag) {
            queue.push_back(Gi);
        }
    }

    for(auto G: queue) {
        if(PQ.size() < k || getScore(G, D) > getScore(PQ.top(), D)) {
            constructGroups(G, l, k, PQ, D);
        } else {
            break;
        }
    }

    while(PQ.size()) {
        topKSkylineGroups.push_back(PQ.top());
        PQ.pop();
    }

    return topKSkylineGroups;
}


Point f_min(Group grp){
    if(grp.points.size()==0){
        cout<<"f_min error"<<endl;
        exit(1);
    }
    vector<int> new_dim;
    for(auto g : grp.points){
        if(new_dim.empty()){
            new_dim = g.dims;
        } else {
            for(int i=0;i<new_dim.size();i++){
                new_dim[i] = min(new_dim[i], g.dims[i]);
            }
        }
    }
    Point Q(new_dim, "Q");
    return Q;
}

vector<Group> skyline_groups_min(vector<Group> sky, vector<Group> temp){
    vector<Group> grp;
    for(auto g1 : temp){
        bool flag = 0;
        Point p1 = f_min(g1);
        for(auto g2 : temp){
            Point p2 = f_min(g2);
            if(p2.dominates(p1)==1){
                flag = 1;
                break;
            }
        }
        if(!flag){
            grp.push_back(g1);
        }
    }
    if(sky.size() == 0){
        return grp;
    }
    vector<Group> skyline;
    for(auto s : sky){
        bool flag = 0;
        for(auto g : grp){
            auto p1 = f_min(s);
            auto p2 = f_min(g);
            if(p2.dominates(p1)==1){
                flag = 1;
                break;
            }
        }
        if(!flag){
            skyline.push_back(s);
        }
    }
    for(auto g : grp){
        bool flag = 0;
        for(auto s : sky){
            auto p1 = f_min(s);
            auto p2 = f_min(g);
            if(p1.dominates(p2)==1){
                flag = 1;
                break;
            }
        }
        if(!flag){
            skyline.push_back(g);
        }
    }
    return skyline;
}

vector<Point> compute_theta_v(Group grp, vector<Point> D){
    Point p = f_min(grp);
    vector<Point> theta;
    for(auto p2 : D){
        if(p2.dominates(p) == 1){
            theta.push_back(p2);
        } else if(p2.dims == p.dims){
            theta.push_back(p2);
        }
    }
    return theta;
}
void generate_l_group(int i,vector<Point> &theta,int l,vector<Group> &grp,Group curgrp){
	if(l==0&&i==theta.size()){
		grp.push_back(curgrp);
		return;
	}
	if(l!=0&&i==theta.size()){
		return;
	}
	generate_l_group(i+1,theta,l,grp,curgrp);
	curgrp.points.insert(theta[i]);
	generate_l_group(i+1,theta,l-1, grp,curgrp);
}

bool compare_group(Group grp1, Group grp2){
    if(grp1.points.size() != grp2.points.size()){
        return 0;
    }
    vector<Point> v1, v2;
    for(auto g : grp1.points) v1.push_back(g);
    for(auto g : grp2.points) v2.push_back(g);

    bool flag = 0;
    for(int i=0;i<v1.size();i++){
        if(v1[i].dims != v2[i].dims){
            flag = 1;
        }
    }
    if(flag) return 0;
    return 1;
}

vector<Group> TKD_MIN(vector<Point> D, int k, int l) {
    D = inputPruning(D, l);
    priority_queue<Group, vector<Group>, comp> PQ;
    vector<Group> Sky[D.size() + 1][l + 1];
    vector<Group> topKSkylineGroups;

    for(int j = 1;j<=D.size();j++){
        for(int i = min(j,l); i>=1;i--){
            if(i==1){
                Group g;
                g.points.insert(D[j-1]);
                vector<Group> dj;
                dj.push_back(g);
                Sky[j][1] = skyline_groups_min(Sky[j-1][1], dj);
            } else {
                vector<Group> temp = Sky[j-1][i-1];
                for(int k=0;k<temp.size();k++){
                    temp[k].points.insert(D[j-1]);
                }
                Sky[j][i] = skyline_groups_min(Sky[j-1][i], temp);
            }
        }
    }
    map<vector<int>,int> mp;
    vector<Group> skyl;
    for(auto g: Sky[D.size()][l]){
        Point p = f_min(g);
        vector<int> dims=p.dims;
        if(mp[dims]==0){
            mp[dims] = 1;
            skyl.push_back(g);
        }
    }
    for(auto g : skyl){
        vector<Point> theta = compute_theta_v(g, D);
        Group theta_g;
        for(auto v : theta) theta_g.points.insert(v);
            if(!PQ.empty()){
                if(PQ.size() < k || getScore(theta_g, D) > getScore(PQ.top(), D)){
                    vector<Group> l_grp;
                    Group curgrp;
                    generate_l_group(0, theta, l, l_grp, curgrp);
                    for(auto g: l_grp){
                        int tm=0;
                        vector<Group> rmvec;
                        
                        while(!PQ.empty()){
                            if(compare_group(g, PQ.top())){
                                tm=1;
                            }
                            rmvec.push_back(PQ.top());
                            PQ.pop();
                        }
                        for(auto x:rmvec){
                            PQ.push(x);
                        }
                        if(PQ.size() < k && tm==0){
                            PQ.push(g);
                        }
                        else if(getScore(g,D)>getScore(PQ.top(),D)&&tm==0){
                            PQ.pop();
                            PQ.push(g);
                        }
                    }
                }
            } else {
                vector<Group> l_grp;
                Group curgrp;
                generate_l_group(0, theta, l, l_grp, curgrp);
                for(auto g: l_grp){
                        int tm=0;
                        vector<Group> rmvec;
                        while(!PQ.empty()){
                            if(compare_group(g, PQ.top())){
                                tm=1;
                            }
                            rmvec.push_back(PQ.top());
                            PQ.pop();
                        }
                        for(auto x:rmvec){
                            PQ.push(x);
                        }
                    if(PQ.size() < k && tm==0){
                        PQ.push(g);
                    }
                    else if(getScore(g,D)>getScore(PQ.top(),D)&&tm==0){
                        PQ.pop();
                        PQ.push(g);
                    }
                }
            }
        
    }
    while(PQ.size()) {
        topKSkylineGroups.push_back(PQ.top()); PQ.pop();
    }
    reverse(topKSkylineGroups.begin(),topKSkylineGroups.end());
    return topKSkylineGroups;
}

int main() {
    // auto topKSkylineGroups = TKD_permutation(D, 4, 3);
    auto topKSkylineGroups = TKD_MAX(D, 4, 2);
    // auto topKSkylineGroups = TKD_SUM(D, 4, 2);
    // auto topKSkylineGroups = TKD_MIN(D, 7, 2);

    for(auto g : topKSkylineGroups) {
        for(auto p : g.points) {
            cout << p.name << ": ";
            for(auto x : p.dims) cout << x << ",";
            cout << " ";
        }
        cout << getScore(g, D) << "\n";
    }
}