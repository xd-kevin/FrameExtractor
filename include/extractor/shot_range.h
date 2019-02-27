//
// Created by shawntang(唐兴) on 2019-02-15.
//

#ifndef FRAMEEXTRACTOR_RANGE_H
#define FRAMEEXTRACTOR_RANGE_H

#include <vector>
using namespace std;

namespace extractor {
    class Range {
    public:
        int start;
        int end;
        std::vector<int> v_idx;
    public:
        Range(int s, int e):start(s),end(e){};

        inline int length() const {return end-start+1;};

        inline void print() const {
            printf("range(%d:%d) (%d) [", start,end,length());
            for(size_t i = 0; i < v_idx.size(); i++){
                printf("%d", v_idx[i]);
                if(i+1 < v_idx.size())
                    printf(",");
            }
            printf("]\n");
        };
    };

    class ShotRange: public Range{
    public:
        vector<Range> v_range;

    public:
        ShotRange(int s, int e): Range(s,e) {};

        inline void print() const{
            printf("shot_range(%d:%d)(%d)[", start, end, length());
            for(size_t i = 0; i < v_idx.size(); i++){
                printf("%d",v_idx[i]);
                if(i+1< v_idx.size())
                    printf(",");
            }
            printf("]\n");

            for(size_t i=0; i < v_range.size(); i++){
                printf(" sub_");
                v_range[i].print();
            }
        }
    };
}


#endif //FRAMEEXTRACTOR_RANGE_H
