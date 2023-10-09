//
// Created by rainbowwing on 2023/8/26.
//

#ifndef THICKEN2_FLAG_PARSER_H
#define THICKEN2_FLAG_PARSER_H
#include <gflags/gflags.h>

DEFINE_int32(m, 1,
             "This arg is a integer means result mode which value can be chose in 1,2,3. \n mode 1 is to remeshing after executing the algorithm\n; mode 2 is cancel the step of remeshing.");

DEFINE_double(l, 1.30,
              "This arg is a double which value indicates how many times the maximum offset distance is the ideal offset distance.");

DEFINE_double(s, 1.0,
              "This arg is a double which value indicates how many times the minimum offset distance is the ideal offset distance.");
DEFINE_int32(t, 12, "thread num please set this value depend the cpu of you device.");
DEFINE_string(f, "", "file name, it must be like *.obj2");
DEFINE_int32(i, 1,
             "This arg is a integer means running mode which value can be chose in 1,2. \n 1 is offsetting to the outside of the mesh; 2 is offsetting to the inside of the mesh.");

DEFINE_double(d, -1,
             "The absolute distance for offset. If you want to perform variable offset, please do not use this parameter, but use the obj2 file.");

DEFINE_double(L, -1,
              "set tetwild argument -l");

DEFINE_double(E, -1,
              "set tetwild argument -e");

//DEFINE_double(e, 1e-4,
//              "This arg is a double means the eps. When the distance of two points is smaller than eps, we will regard these two point as coinciding. ");
int result_mode;
string input_filename;
double tetwild_l = -1;
double tetwild_e = -1;
void flag_parser() {
    result_mode = FLAGS_m;
    cout << "result_mode is " << result_mode << endl;
    running_mode = FLAGS_i;
    if (running_mode != 1 && running_mode != 2) {
        cout << "running_mode must equal to 1 or 2" << endl;
        exit(0);
    }
    min_distance_limit = FLAGS_s;
    max_distance_limit = FLAGS_l;
    thread_num = FLAGS_t;
    input_filename = FLAGS_f;
    myeps = 0;
    tetwild_l = FLAGS_L;
    tetwild_e = FLAGS_E;
    default_move = FLAGS_d;

}


#endif //THICKEN2_FLAG_PARSER_H
