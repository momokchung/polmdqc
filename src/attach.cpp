////////////////////////////////////////////////////
//                                                //
//  attach.cpp  --  setup of connectivity arrays  //
//                                                //
////////////////////////////////////////////////////

// "attach" generates lists of 1-3, 1-4 and 1-5 connectivities
// starting from the previously determined list of attached
// atoms (ie, 1-2 connectivity)


#include "atoms.h"
#include "attach.h"
#include "couple.h"
#include "fatal.h"
#include "sort.h"

void attach()
{
    // perform dynamic allocation of some global arrays
    int maxn13 = 3 * maxval;
    int maxn14 = 9 * maxval;
    int maxn15 = 27 * maxval;
    n13.resize(n, 0);
    n14.resize(n, 0);
    n15.resize(n, 0);
    i13.resize(n, std::vector<int>(maxn13));
    i14.resize(n, std::vector<int>(maxn14));
    i15.resize(n, std::vector<int>(maxn15));

    // loop over all atoms finding all the 1-3 relationships;
    // note "n12" and "i12" have already been setup elsewhere
    for (int i = 0; i < n; i++) {
        std::vector<int> i13tmp;
        for (int j = 0; j < n12[i]; j++) {
            int jj = i12[i][j];
            for (int k = 0; k < n12[jj]; k++) {
                int kk = i12[jj][k];
                if (kk == i) goto label_10;
                for (int m = 0; m < n12[i]; m++){
                    if (kk == i12[i][m])  goto label_10;
                }
                i13tmp.push_back(kk);
                label_10:
                continue;
            }
        }
        int n13tmp = i13tmp.size();
        if (n13tmp > maxn13) {
            printf("\n ATTACH  --  Too many 1-3 Connected Atoms Attached to Atom%d", i);
            fatal();
        }
        sort(i13tmp);
        // std::sort(i13[i], i13[i] + n13[i]);
    }
}


// c

// c
// c     loop over all atoms finding all the 1-4 relationships
// c
//       do i = 1, n
//          n14(i) = 0
//          do j = 1, n13(i)
//             jj = i13(j,i)
//             do k = 1, n12(jj)
//                kk = i12(k,jj)
//                if (kk .eq. i)  goto 30
//                do m = 1, n12(i)
//                   if (kk .eq. i12(m,i))  goto 30
//                end do
//                do m = 1, n13(i)
//                   if (kk .eq. i13(m,i))  goto 30
//                end do
//                n14(i) = n14(i) + 1
//                i14(n14(i),i) = kk
//    30          continue
//             end do
//          end do
//          if (n14(i) .gt. maxn14) then
//             write (iout,40)  i
//    40       format (/,' ATTACH  --  Too many 1-4 Connected Atoms',
//      &                 ' Attached to Atom',i6)
//             call fatal
//          end if
//          call sort8 (n14(i),i14(1,i))
//       end do
// c
// c     loop over all atoms finding all the 1-5 relationships
// c
//       do i = 1, n
//          n15(i) = 0
//          do j = 1, n14(i)
//             jj = i14(j,i)
//             do k = 1, n12(jj)
//                kk = i12(k,jj)
//                if (kk .eq. i)  goto 50
//                do m = 1, n12(i)
//                   if (kk .eq. i12(m,i))  goto 50
//                end do
//                do m = 1, n13(i)
//                   if (kk .eq. i13(m,i))  goto 50
//                end do
//                do m = 1, n14(i)
//                   if (kk .eq. i14(m,i))  goto 50
//                end do
//                n15(i) = n15(i) + 1
//                i15(n15(i),i) = kk
//    50          continue
//             end do
//          end do
//          if (n15(i) .gt. maxn15) then
//             write (iout,60)  i
//    60       format (/,' ATTACH  --  Too many 1-5 Connected Atoms',
//      &                 ' Attached to Atom',i6)
//             call fatal
//          end if
//          call sort8 (n15(i),i15(1,i))
//       end do
//       return
//       end
