//
// Created by shawntang(唐兴) on 2019-02-14.
//

#ifndef FRAMEEXTRACTOR_EXTRACTOR_H
#define FRAMEEXTRACTOR_EXTRACTOR_H

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <getopt.h>
#include <sys/stat.h>
#include <algorithm>

#include <string>


using namespace std;

//program options

struct extractor_params{
    string in_video;
    string out_dir;

    int step_sz; //frame subsampling
    double fltr_begin_sec; //filter begin
    double fltr_end_sec; //filter end
    double invalid_wnd; //window for dropping neighbor frames of low-quality ones
    bool sample; //if enabled, random sample

    extractor_params():
    out_dir("./output"),
    step_sz(1),
    fltr_begin_sec(-1.0),
    fltr_end_sec(-1.0),
    invalid_wnd(0.15),
    sample(false)
    {};
};

inline void extractor_parse_params(int argc, char** argv, extractor_params& opt)
{
    static struct option long_options[] = {
            {"in_video",required_argument,0,'i'},
            {"out_dir", required_argument,0,'o'},
            {"step", required_argument, 0, 's'},
            {"fltr_begin_sec", required_argument, 0, 'b'},
            {"fltr_end_sec", required_argument, 0, 'e'},
            {"invalid_wnd", required_argument, 0, 'w'},
            {"sample", no_argument, 0, 'S'},
            {0,0,0,0}
    };

    while (true){
        int opt_idx = 0;
        int c = getopt_long(argc, argv,"i:o:s:b:e:w:"
                                       "S", long_options, &opt_idx);
        if(c == -1) break;
        switch (c) {
            case 'i': opt.in_video = optarg; break;
            case 'o': opt.out_dir = optarg; break;
            case 's': opt.step_sz = atoi(optarg); break;
            case 'b': opt.fltr_begin_sec = atof(optarg); break;
            case 'e': opt.fltr_end_sec = atof(optarg); break;
            case 'w': opt.invalid_wnd = atof(optarg); break;
            case 'S': opt.sample = true; break;
        }
    }
    mkdir( opt.out_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}
inline void usage()
{
    extractor_params opt;
    printf("USAGE: hecate -i infile [options]\n");
    printf("\n");
    printf("  -i  --in_video      (string)    Input video file\n");
    printf("  -o  --out_dir       (string)    Output directory (%s)\n", opt.out_dir.c_str());
    printf("  -S  --sample          (bool,option)       Frame subsampling\n");
    exit(-1);
}

void run_extractor(extractor_params& opt, vector<int>& v_frm_idx);




#endif //FRAMEEXTRACTOR_EXTRACTOR_H
