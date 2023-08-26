//
// Created by rainbowwing on 2023/8/26.
//

#ifndef THICKEN2_FLAG_PARSER_H
#define THICKEN2_FLAG_PARSER_H
#include <gflags/gflags.h>

DEFINE_int32(m, 3,
             "This arg is a integer means result mode which value can be chose in 1,2,3. \n mode 1  without step 6;mode 2 is with building mesh;mode 3 is with building mesh and remeshing.");
DEFINE_double(d, 0.5,
              "This arg is a double means the length running invariable thickening limited in 0.1~1.5. Example -d 0.5, Indicating that the offset distance is 0.75 times the average mesh edge length.");
DEFINE_double(l, 2.0,
              "This arg is a double which value indicates how many times the maximum offset distance is the ideal offset distance limited in 1.5~2.7. .You can set it is 2.0");
DEFINE_double(g, -1,
              "This arg is a double means the length of edge length of each cell in grid. If you can't calculate a length with better performance, it can be passed. Then it will use the default value.");
DEFINE_int32(t, 12, "thread num please set this value depend the cpu of you device.");
DEFINE_string(f, "", "file name, which can choose *.obj2 or *.obj.");
DEFINE_int32(i, 1,
             "This arg is a integer means running mode which value can be chose in 1,2. \n 1 is offsetting to the outside of the mesh; 2 is offsetting to the inside of the mesh.");
DEFINE_double(e, 1e-5,
              "This arg is a double means the eps. When the distance of two points is smaller than eps, we will regard these two point as coinciding. ");
DEFINE_bool(s, false,
            "This arg is a bool value which means the program will skip some cell which is most likely useless. It can improve performance, but may cause holes in the result. We suggest not use this function.");
void flag_parser() {
    int result_mode = FLAGS_m;
    cout << "result_mode is " << result_mode << endl;
    int running_mode = FLAGS_i;
    if (running_mode != 1 && running_mode != 2) {
        cout << "running_mode must equal to 1 or 2" << endl;
        exit(0);
    }
}


#endif //THICKEN2_FLAG_PARSER_H
