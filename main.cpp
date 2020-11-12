#include<bits/stdc++.h>
using namespace std;


class Point {
public:
    string name;
    vector<int> dims;
    Point(vector<int> dims, string name) : name(name){
        this->dims = dims;
    }
    
    int dominates(Point& other) {
        int a = 0, b = 0, c = 0;
        int d = other.dims.size();
        for(int i = 0; i < d; i++) {
            if(dims[i] < other.dims[i]) b++;
            else if(dims[i] > other.dims[i]) a++; 
            else c++;
        }
        if((a + c) == d) return 1;
        else if((a + c) < d && (b + c) == d) return -1;
        return 0;
    }

    bool operator<(const Point& other) const {
        return 1;
    }
};

vector<Point> D = {Point({10, 0, 0}, "Q1"), Point({1,7,1}, "Q2"), Point({2,3,3}, "Q3"), Point({7,0,0}, "Q4"), Point({1,5,0}, "Q5"), Point({0,0,4}, "Q7"), Point({0,4,1}, "Q6"), Point({1,1,2}, "Q8"), Point({0,1,3}, "Q9"), Point({1,0,3}, "Q10")};

class Group {
public:
    set<Point> points;
    int score;

    Group(set<Point> points) : points(points) {
        score = 0;
    }

};

int getScore(Group G, vector<Point>& D) {
    map<Point, bool> hmap;
    int score = 0;
    for(auto q : G.points) {
        for(auto p : D) {
            if(q.dominates(p) == 1 && !hmap[p]) {
                hmap[p] = true;
                score++;
            }
        }
    }
    return score;
}


struct comp{ 
    bool operator()(const Group& curr, const Group& other) {
        return getScore(curr, D) > getScore(other, D);
    }
};
vector<Point> inputPruning(vector<Point> D) {

} 

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

vector<Point> findResidual(vector<Point>& D, vector<Point>& skyline) {
    map<Point, bool> hmap;
    for(auto s : skyline) hmap[s] = true;

    vector<Point> residual;
    for(auto p : D) {
        if(!hmap[p]) {
            residual.push_back(p);
        }
    }
    return residual;
}


void generateGroups(int pos, int l, int k, priority_queue<Group, vector<Group>, comp>& PQ, vector<Point> skyline, Group candidate) {
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

void revive(int pos, int l, int k, priority_queue<Group, vector<Group>, comp>& PQ, vector<Point> residual, Group candidate) {
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
                if(res == 0 || res == -1) {
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
    auto skyline = findSkyline(D);
    auto residual = findResidual(D, skyline);
    set<Group, comp> candidate, bucket;
    priority_queue<Group, vector<Group>, comp> PQ;
    vector<Group> topKSkylineGroups;

    generateGroups(0, l, k, PQ, skyline, Group(set<Point>()));

    while(PQ.size()) {
        auto G = PQ.top(); PQ.pop();
        // for(auto p : G.points) cout << p.name << " ";
        // cout << "\n";
        if(G.points.size() < l) {
            candidate.insert(G);
        } else {
            bucket.insert(G);
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

int main() {
    auto topKSkylineGroups = TKD_permutation(D, 4, 3);

    for(auto g : topKSkylineGroups) {
        for(auto p : g.points) {
            cout << p.name << ": ";
            for(auto x : p.dims) cout << x << ",";
            cout << " ";
        }
        cout << getScore(g, D) << "\n";
    }
}