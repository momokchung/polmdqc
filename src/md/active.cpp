// Author: Moses KJ Chung
// Year:   2023

#include "active.h"
#include "atoms.h"
#include "gettext.h"
#include "inform.h"
#include "keys.h"
#include "upcase.h"
#include "usage.h"
#include <sstream>
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////
//                                            //
//  active  --  set the list of active atoms  //
//                                            //
////////////////////////////////////////////////

// "active" sets the list of atoms that are used during
// coordinate manipulation or potential energy calculations

void active()
{
    int next;
    int nmobile,nfixed;
    int center,nsphere;
    real xr,yr,zr;
    real xcenter,ycenter,zcenter;
    real radius,radius2,dist2;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;
    bool header;

    // allocation and initialization of some global arrays
    iuse.allocate(n);
    use.allocate(n+1);

    // allocation and initialization of some local arrays
    std::vector<int> mobile;
    std::vector<int> fixed;
    mobile.resize(n);
    fixed.resize(n);

    // set defaults for the numbers and lists of active atoms
    nuse = n;
    use[0] = false;
    for (int i = 0; i < n; i++) {
        iuse[i] = -1;
        use[i+1] = true;
    }
    nmobile = 0;
    nfixed = 0;
    for (int i = 0; i < n; i++) {
        mobile[i] = 0;
        fixed[i] = 0;
    }
    nsphere = 0;

    // get any keywords containing active atom parameters
    for (int j = 0; j < nkey; j++) {
        record = keyline[j];
        next = 0;
        gettext(record, keyword, next);
        upcase(keyword);
        string = record.substr(next);
        iss.clear();
        iss.str(string);

        // get any lists of atoms whose coordinates are active
        if (keyword == "ACTIVE") {
            int mobileInt = 0;
            int counter = 0;
            while (iss >> mobileInt) {
                mobile[nmobile + counter] = mobileInt;
                counter++;
            }
            while (mobile[nmobile] != 0) {
                nmobile++;
            }
        }

        // get any lists of atoms whose coordinates are inactive
        else if (keyword == "INACTIVE") {
            int fixedInt = 0;
            int counter = 0;
            while (iss >> fixedInt) {
                fixed[nfixed + counter] = fixedInt;
                counter++;
            }
            while (fixed[nfixed] != 0) {
                nfixed++;
            }
        }

        // get the center and radius of the sphere of active atoms
        else if (keyword == "ACTIVE-SPHERE") {
            center = -1;
            xcenter = 0.;
            ycenter = 0.;
            zcenter = 0.;
            radius = 0.;
            iss >> xcenter >> ycenter >> zcenter >> radius;
            if (radius == 0.) {
                iss.clear();
                iss.str(string);
                if (!(iss >> center >> radius)) goto label_60;
                xcenter = x[center-1];
                ycenter = y[center-1];
                zcenter = z[center-1];
            }
            nsphere++;
            if (nsphere == 1) {
                nuse = 0;
                for (int i = 0; i < n; i++) {
                    use[i+1] = false;
                }
                if (verbose) {
                    std::string blank6(6, ' ');
                    std::string blank11(11, ' ');
                    std::string blank12(12, ' ');
                    printf("\n Spheres used to Select Active Atoms :\n");
                    printf("\n   Atom Center%sCoordinates%sRadius%s# Active Atoms\n", blank11.c_str(), blank12.c_str(), blank6.c_str());
                }
            }
            radius2 = radius * radius;
            for (int i = 0; i < n; i++) {
                if (!use[i+1]) {
                    xr = x[i] - xcenter;
                    yr = y[i] - ycenter;
                    zr = z[i] - zcenter;
                    dist2 = xr*xr + yr*yr + zr*zr;
                    if (dist2 <= radius2) {
                        nuse++;
                        use[i+1] = true;
                    }
                }
            }
            if (verbose) {
                printf("  %8d      %9.2f%9.2f%9.2f  %9.2f       %8d\n", center, xcenter, ycenter, zcenter, radius, nuse);
            }
            label_60:
            continue;
        }
    }

    // remove active or inactive atoms not in the system
    header = true;
    for (int i = 0; i < n; i++) {
        if (std::abs(mobile[i]) > n) {
            mobile[i] = 0;
            if (header) {
                header = false;
                printf("\n ACTIVE  --  Warning, Illegal Atom Number in ACTIVE Atom List\n");
            }
        }
    }
    header = true;
    for (int i = 0; i < n; i++) {
        if (std::abs(fixed[i]) > n) {
            fixed[i] = 0;
            if (header) {
                header = false;
                printf("\n ACTIVE  --  Warning, Illegal Atom Number in ACTIVE Atom List\n");
            }
        }
    }

    // set active atoms to those marked as not inactive
    int i = 0;
    while (fixed[i] != 0) {
        if (fixed[i] > 0) {
            int j = fixed[i]-1;
            if (use[j+1]) {
                use[j+1] = false;
                nuse--;
            }
            i++;
        }
        else {
            for (int j = abs(fixed[i])-1; j < abs(fixed[i+1]); j++) {
                if (use[j+1]) {
                    use[j+1] = false;
                    nuse--;
                }
            }
            i += 2;
        }
    }

    // set active atoms to only those marked as active
    i = 0;
    while (mobile[i] != 0) {
        if (i == 0) {
            nuse = 0;
            for (int j = 0; j < n; j++) {
                use[j+1] = false;
            }
        }
        if (mobile[i] > 0) {
            int j = mobile[i]-1;
            if (!use[j+1]) {
                use[j+1] = true;
                nuse += 1;
            }
            i += 1;
        }
        else {
            for (int j = abs(mobile[i])-1; j < abs(mobile[i+1]); j++) {
                if (!use[j+1]) {
                    use[j+1] = true;
                    nuse += 1;
                }
            }
            i += 2;
        }
    }

    // use logical array to set the list of active atoms
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (use[i+1]) {
            iuse[j] = i;
            j += 1;
        }
    }

    // output the final list of the active atoms
    if (debug and nuse>0 and nuse<n) {
        printf("\n List of Active Atoms for Energy Calculations :\n\n   ");
        for (int i = 0; i < nuse; i++) {
            printf("%7d", iuse[i]+1);
            if ((i + 1) % 10 == 0) printf("\n   ");
        }
        printf("\n");
    }
}
}
