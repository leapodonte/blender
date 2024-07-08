#ifndef ORTHOLAB_MODEL_H
#define ORTHOLAB_MODEL_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdint.h>
#include "bindings_types.h"

    OrthoLabModel_t *new_ortholab();
    void delete_ortholab(OrthoLabModel_t *point_list);






    // test cases
    int ortholab_model_test_db();   
    int ortholab_model_test_prj();   

#ifdef __cplusplus
}
#endif

#endif
