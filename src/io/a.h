// Part of SimCardiac Project.

#pragma once
#include <petscdm.h>
#include <petscksp.h>
#include <petscdmstag.h> /* Includes petscdmproduct.h */


static PetscErrorCode PrintReferenceSolution_v1(
    DM dmSol, Vec *pSolRef)
{
  PetscInt       startx, starty, startgx, startgy, nx, ny, nExtra[2], ex, ey, lx, ly, gx, gy;
  PetscInt       iuy, iux, ip, iprev, icenter, ipnext;
  PetscScalar ***arrSol, **cArrX, **cArrY;
  Vec            solRefLocal;

  PetscFunctionBeginUser;
  PetscCall(DMCreateGlobalVector(dmSol, pSolRef));
  PetscCall(DMGetLocalVector(dmSol, &solRefLocal));

  /* Obtain indices to use in the raw arrays */
  PetscCall(DMStagGetLocationSlot(dmSol, DMSTAG_DOWN, 0, &iuy)); /* iuy 表示单元下边界的第 0 个分量 */
  PetscCall(DMStagGetLocationSlot(dmSol, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dmSol, DMSTAG_ELEMENT, 0, &ip));

  /* Use high-level convenience functions to get raw arrays and indices for 1d coordinates */
  PetscCall(DMStagGetProductCoordinateArraysRead(dmSol, &cArrX, &cArrY, NULL)); /* cArrX[ex][icenter] 表示坐标 */
  PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, DMSTAG_LEFT, &iprev));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, DMSTAG_RIGHT, &ipnext));

  PetscCall(DMStagVecGetArray(dmSol, solRefLocal, &arrSol));
  PetscCall(DMStagGetCorners(dmSol, &startx, &starty, NULL, &nx, &ny, NULL, &nExtra[0], &nExtra[1], NULL));
  PetscCall(DMStagGetGhostCorners(dmSol, &startgx, &startgy, NULL, &gx, &gy, NULL));
  PetscCall(DMStagGetLocalSizes(dmSol, &lx, &ly, NULL));
  printf("%p\n", (void*)(*arrSol));
  printf("startx: %d, starty: %d, nx: %d, ny: %d,\n"
         "startgx: %d, startgy: %d, gx: %d, gy: %d,\n"
         "    lx: %d, ly: %d, extra_x: %d, extra_y: %d\n", startx, starty, nx, ny, startgx, startgy, gx, gy, lx, ly, nExtra[0], nExtra[1]);
//   for (ey = starty; ey < starty + ny + nExtra[1]; ++ey) {
//     for (ex = startx; ex < startx + nx + nExtra[0]; ++ex) {
//         printf("Element (%d, %d): y : (%f, %f), x : (%f, %f), p : (%f, %f)\n", ex, ey, 
//             cArrX[ex][icenter], cArrY[ey][iprev],
//             cArrX[ex][iprev], cArrY[ey][icenter],
//             cArrX[ex][icenter], cArrY[ey][ipnext]);
//       if (!(ey < starty + ny && ex < startx + nx)) { /* Don't fill on the dummy elements (though you could, and these values would just be ignored) */
//         printf("Dummy Element (%d, %d)\n", ex, ey);
//         arrSol[ey][ex][ip] = 0.0;
//         printf("  p = %g\n", arrSol[ey][ex][ip]);
//       }
//     }
//   }
  PetscCall(DMStagVecRestoreArray(dmSol, solRefLocal, &arrSol));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dmSol, &cArrX, &cArrY, NULL));
  PetscCall(DMLocalToGlobal(dmSol, solRefLocal, INSERT_VALUES, *pSolRef));
  PetscCall(DMRestoreLocalVector(dmSol, &solRefLocal));
  PetscFunctionReturn(PETSC_SUCCESS);
}

