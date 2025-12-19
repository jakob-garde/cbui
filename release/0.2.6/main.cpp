// #define ENABLE_GLFW

#include "jg_cbui.h"


int main (int argc, char **argv) {
    printf("Check cbui %u.%u.%u build and version", CBUI_VERSION_MAJOR, CBUI_VERSION_MINOR, CBUI_VERSION_PATCH);

    BaselayerAssertVersion(0, 2, 6);
    CbuiAssertVersion(0, 2, 6);

    printf(" - OK\n");
}
