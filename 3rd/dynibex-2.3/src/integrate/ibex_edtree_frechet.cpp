/* ============================================================================
 * D Y N I B E X - Definition of the Trees based on Frechet derivatives
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Author(s) of modification   : Julien Alexandre dit Sandretto and Olivier Mullier
 * Modified    : Jul 25, 2018
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#include "ibex_Affine2.h" 
#include "ibex_Affine2Vector.h"
#include "ibex_Function.h"
#include "ibex_edtree_frechet.h"

#include <vector>
namespace ibex{

Affine2 edtree_frechet::get_derivatives(int order, Affine2Vector y, int j){
    int j1,j2,j3,j4,j5;
    switch (order){
      case 0 :
	return y[j];
case 1:
{
  Affine2 res1(0.0);
 std::vector<int> key0;
  // Tree computation

  key0.clear();
  key0.push_back (j);
  res1 += edfr->eval_frechet(0, key0, y);

  return res1;
}

case 2:
{
  Affine2 res1(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp1 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp2 = edfr->eval_frechet(1, key0, y);

    // Tree computation
    res1 += temp1 * temp2;
  }

  return res1;
}

case 3:
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp2 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp5 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp1 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp3 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp4 = edfr->eval_frechet(1, key0, y);

      // Tree computation
      res1 += temp1 * temp2 * temp3;
      res2 += temp1 * temp4 * temp5;
    }
  }

  return res2 + res1;
}

case 4:
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
 std::vector<int> key3;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp3 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp8 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp2 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp6 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp9 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp1 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp4 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp7 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp5 = edfr->eval_frechet(1, key0, y);

        // Tree computation
        res1 += temp1 * temp2 * temp3 * temp4;
        res2 += temp1 * temp5 * temp3 * temp6;
        res3 += temp1 * temp2 * temp7 * temp8;
        res4 += temp1 * temp5 * temp9 * temp8;
      }

    }


  }

  return 3 * (res2) + res4 + res3 + res1;
}

case 5:
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  Affine2 res5(0.0);
  Affine2 res6(0.0);
  Affine2 res7(0.0);
  Affine2 res8(0.0);
  Affine2 res9(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
 std::vector<int> key3;
 std::vector<int> key4;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp4 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp14 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp3 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp9 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp16 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp2 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp7 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp15 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j3);
        Affine2 temp12 = edfr->eval_frechet(1, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp10 = edfr->eval_frechet(1, key0, y);

        for (int j4 = 0; j4 < nbvar; ++j4) {
          // Common Factor

          key0.clear();
          key0.push_back (j4);
          Affine2 temp1 = edfr->eval_frechet(0, key0, y);

          key0.clear();
          key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp5 = edfr->eval_frechet(4, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp13 = edfr->eval_frechet(3, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp8 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j4);
          Affine2 temp11 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j3); key0.push_back (j4);
          Affine2 temp6 = edfr->eval_frechet(1, key0, y);

          // Tree computation
          res1 += temp1 * temp2 * temp3 * temp4 * temp5;
          res2 += temp1 * temp6 * temp3 * temp4 * temp7;
          res3 += temp1 * temp2 * temp8 * temp4 * temp9;
          res4 += temp1 * temp6 * temp10 * temp4 * temp9;
          res5 += temp1 * temp11 * temp2 * temp12 * temp9;
          res6 += temp1 * temp2 * temp3 * temp13 * temp14;
          res7 += temp1 * temp6 * temp3 * temp15 * temp14;
          res8 += temp1 * temp2 * temp8 * temp16 * temp14;
          res9 += temp1 * temp6 * temp10 * temp16 * temp14;
        }
      }
    }
  }

  return 6 * (res2) + 4 * (res4 + res3) + 3 * (res7 + res5) + res9 + res8 + res6 + res1;
}
case 6:
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  Affine2 res5(0.0);
  Affine2 res6(0.0);
  Affine2 res7(0.0);
  Affine2 res8(0.0);
  Affine2 res9(0.0);
  Affine2 res10(0.0);
  Affine2 res11(0.0);
  Affine2 res12(0.0);
  Affine2 res13(0.0);
  Affine2 res14(0.0);
  Affine2 res15(0.0);
  Affine2 res16(0.0);
  Affine2 res17(0.0);
  Affine2 res18(0.0);
  Affine2 res19(0.0);
  Affine2 res20(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
 std::vector<int> key3;
 std::vector<int> key4;
 std::vector<int> key5;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp5 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp21 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp4 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp15 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp24 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp3 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp10 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp23 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j3);
        Affine2 temp19 = edfr->eval_frechet(1, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp17 = edfr->eval_frechet(1, key0, y);

        for (int j4 = 0; j4 < nbvar; ++j4) {
          // Common Factor

          key0.clear();
          key0.push_back (j4);
          Affine2 temp2 = edfr->eval_frechet(0, key0, y);

          key0.clear();
          key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp8 = edfr->eval_frechet(4, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp22 = edfr->eval_frechet(3, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp16 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j4);
          Affine2 temp13 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j3); key0.push_back (j4);
          Affine2 temp11 = edfr->eval_frechet(1, key0, y);

          for (int j5 = 0; j5 < nbvar; ++j5) {
            // Common Factor

            key0.clear();
            key0.push_back (j5);
            Affine2 temp1 = edfr->eval_frechet(0, key0, y);

            key0.clear();
            key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp6 = edfr->eval_frechet(5, key0, y);

            key0.clear();
            key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp20 = edfr->eval_frechet(4, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp14 = edfr->eval_frechet(3, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp18 = edfr->eval_frechet(2, key0, y);

            key0.clear();
            key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp9 = edfr->eval_frechet(2, key0, y);

            key0.clear();
            key0.push_back (j3); key0.push_back (j5);
            Affine2 temp12 = edfr->eval_frechet(1, key0, y);

            key0.clear();
            key0.push_back (j4); key0.push_back (j5);
            Affine2 temp7 = edfr->eval_frechet(1, key0, y);

            // Tree computation
            res1 += temp1 * temp2 * temp3 * temp4 * temp5 * temp6;
            res2 += temp1 * temp7 * temp3 * temp4 * temp5 * temp8;
            res3 += temp1 * temp2 * temp9 * temp4 * temp5 * temp10;
            res4 += temp1 * temp7 * temp11 * temp4 * temp5 * temp10;
            res5 += temp1 * temp12 * temp2 * temp13 * temp5 * temp10;
            res6 += temp1 * temp2 * temp3 * temp14 * temp5 * temp15;
            res7 += temp1 * temp7 * temp3 * temp16 * temp5 * temp15;
            res8 += temp1 * temp2 * temp9 * temp17 * temp5 * temp15;
            res9 += temp1 * temp7 * temp11 * temp17 * temp5 * temp15;
            res10 += temp1 * temp2 * temp18 * temp3 * temp19 * temp15;
            res11 += temp1 * temp7 * temp13 * temp3 * temp19 * temp15;
            res12 += temp1 * temp2 * temp3 * temp4 * temp20 * temp21;
            res13 += temp1 * temp7 * temp3 * temp4 * temp22 * temp21;
            res14 += temp1 * temp2 * temp9 * temp4 * temp23 * temp21;
            res15 += temp1 * temp7 * temp11 * temp4 * temp23 * temp21;
            res16 += temp1 * temp12 * temp2 * temp13 * temp23 * temp21;
            res17 += temp1 * temp2 * temp3 * temp14 * temp24 * temp21;
            res18 += temp1 * temp7 * temp3 * temp16 * temp24 * temp21;
            res19 += temp1 * temp2 * temp9 * temp17 * temp24 * temp21;
            res20 += temp1 * temp7 * temp11 * temp17 * temp24 * temp21;
          }
        }
      }
    }
  }

  return 15 * (res7 + res5) + 10 * (res11 + res10 + res4 + res3 + res2) + 6 * (res13) + 5 * (res9 + res8 + res6) + 4 * (res15 + res14) + 3 * (res18 + res16) + res20 + res19 + res17 + res12 + res1;
}

case 7:
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  Affine2 res5(0.0);
  Affine2 res6(0.0);
  Affine2 res7(0.0);
  Affine2 res8(0.0);
  Affine2 res9(0.0);
  Affine2 res10(0.0);
  Affine2 res11(0.0);
  Affine2 res12(0.0);
  Affine2 res13(0.0);
  Affine2 res14(0.0);
  Affine2 res15(0.0);
  Affine2 res16(0.0);
  Affine2 res17(0.0);
  Affine2 res18(0.0);
  Affine2 res19(0.0);
  Affine2 res20(0.0);
  Affine2 res21(0.0);
  Affine2 res22(0.0);
  Affine2 res23(0.0);
  Affine2 res24(0.0);
  Affine2 res25(0.0);
  Affine2 res26(0.0);
  Affine2 res27(0.0);
  Affine2 res28(0.0);
  Affine2 res29(0.0);
  Affine2 res30(0.0);
  Affine2 res31(0.0);
  Affine2 res32(0.0);
  Affine2 res33(0.0);
  Affine2 res34(0.0);
  Affine2 res35(0.0);
  Affine2 res36(0.0);
  Affine2 res37(0.0);
  Affine2 res38(0.0);
  Affine2 res39(0.0);
  Affine2 res40(0.0);
  Affine2 res41(0.0);
  Affine2 res42(0.0);
  Affine2 res43(0.0);
  Affine2 res44(0.0);
  Affine2 res45(0.0);
  Affine2 res46(0.0);
  Affine2 res47(0.0);
  Affine2 res48(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
 std::vector<int> key3;
 std::vector<int> key4;
 std::vector<int> key5;
 std::vector<int> key6;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp6 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp35 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp5 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp22 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp39 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp4 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp16 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp38 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j3);
        Affine2 temp30 = edfr->eval_frechet(1, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp25 = edfr->eval_frechet(1, key0, y);

        for (int j4 = 0; j4 < nbvar; ++j4) {
          // Common Factor

          key0.clear();
          key0.push_back (j4);
          Affine2 temp3 = edfr->eval_frechet(0, key0, y);

          key0.clear();
          key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp11 = edfr->eval_frechet(4, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp37 = edfr->eval_frechet(3, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp33 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j4);
          Affine2 temp28 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp24 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j4);
          Affine2 temp20 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j3); key0.push_back (j4);
          Affine2 temp18 = edfr->eval_frechet(1, key0, y);

          for (int j5 = 0; j5 < nbvar; ++j5) {
            // Common Factor

            key0.clear();
            key0.push_back (j5);
            Affine2 temp2 = edfr->eval_frechet(0, key0, y);

            key0.clear();
            key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp9 = edfr->eval_frechet(5, key0, y);

            key0.clear();
            key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp36 = edfr->eval_frechet(4, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp23 = edfr->eval_frechet(3, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp31 = edfr->eval_frechet(2, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j5);
            Affine2 temp27 = edfr->eval_frechet(1, key0, y);

            key0.clear();
            key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp17 = edfr->eval_frechet(2, key0, y);

            key0.clear();
            key0.push_back (j3); key0.push_back (j5);
            Affine2 temp14 = edfr->eval_frechet(1, key0, y);

            key0.clear();
            key0.push_back (j4); key0.push_back (j5);
            Affine2 temp12 = edfr->eval_frechet(1, key0, y);

            for (int j6 = 0; j6 < nbvar; ++j6) {
              // Common Factor

              key0.clear();
              key0.push_back (j6);
              Affine2 temp1 = edfr->eval_frechet(0, key0, y);

              key0.clear();
              key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp7 = edfr->eval_frechet(6, key0, y);

              key0.clear();
              key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp34 = edfr->eval_frechet(5, key0, y);

              key0.clear();
              key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp21 = edfr->eval_frechet(4, key0, y);

              key0.clear();
              key0.push_back (j2); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp29 = edfr->eval_frechet(3, key0, y);

              key0.clear();
              key0.push_back (j2); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp32 = edfr->eval_frechet(2, key0, y);

              key0.clear();
              key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp15 = edfr->eval_frechet(3, key0, y);

              key0.clear();
              key0.push_back (j3); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp19 = edfr->eval_frechet(2, key0, y);

              key0.clear();
              key0.push_back (j3); key0.push_back (j6);
              Affine2 temp26 = edfr->eval_frechet(1, key0, y);

              key0.clear();
              key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp10 = edfr->eval_frechet(2, key0, y);

              key0.clear();
              key0.push_back (j4); key0.push_back (j6);
              Affine2 temp13 = edfr->eval_frechet(1, key0, y);

              key0.clear();
              key0.push_back (j5); key0.push_back (j6);
              Affine2 temp8 = edfr->eval_frechet(1, key0, y);

              // Tree computation
              res1 += temp1 * temp2 * temp3 * temp4 * temp5 * temp6 * temp7;
              res2 += temp1 * temp8 * temp3 * temp4 * temp5 * temp6 * temp9;
              res3 += temp1 * temp2 * temp10 * temp4 * temp5 * temp6 * temp11;
              res4 += temp1 * temp8 * temp12 * temp4 * temp5 * temp6 * temp11;
              res5 += temp1 * temp13 * temp2 * temp14 * temp5 * temp6 * temp11;
              res6 += temp1 * temp2 * temp3 * temp15 * temp5 * temp6 * temp16;
              res7 += temp1 * temp8 * temp3 * temp17 * temp5 * temp6 * temp16;
              res8 += temp1 * temp2 * temp10 * temp18 * temp5 * temp6 * temp16;
              res9 += temp1 * temp8 * temp12 * temp18 * temp5 * temp6 * temp16;
              res10 += temp1 * temp2 * temp19 * temp3 * temp20 * temp6 * temp16;
              res11 += temp1 * temp8 * temp14 * temp3 * temp20 * temp6 * temp16;
              res12 += temp1 * temp2 * temp3 * temp4 * temp21 * temp6 * temp22;
              res13 += temp1 * temp8 * temp3 * temp4 * temp23 * temp6 * temp22;
              res14 += temp1 * temp2 * temp10 * temp4 * temp24 * temp6 * temp22;
              res15 += temp1 * temp8 * temp12 * temp4 * temp24 * temp6 * temp22;
              res16 += temp1 * temp13 * temp2 * temp14 * temp24 * temp6 * temp22;
              res17 += temp1 * temp2 * temp3 * temp15 * temp25 * temp6 * temp22;
              res18 += temp1 * temp8 * temp3 * temp17 * temp25 * temp6 * temp22;
              res19 += temp1 * temp2 * temp10 * temp18 * temp25 * temp6 * temp22;
              res20 += temp1 * temp8 * temp12 * temp18 * temp25 * temp6 * temp22;
              res21 += temp1 * temp26 * temp2 * temp27 * temp3 * temp28 * temp16;
              res22 += temp1 * temp2 * temp3 * temp29 * temp4 * temp30 * temp22;
              res23 += temp1 * temp8 * temp3 * temp31 * temp4 * temp30 * temp22;
              res24 += temp1 * temp2 * temp10 * temp20 * temp4 * temp30 * temp22;
              res25 += temp1 * temp8 * temp12 * temp20 * temp4 * temp30 * temp22;
              res26 += temp1 * temp2 * temp32 * temp3 * temp4 * temp33 * temp22;
              res27 += temp1 * temp8 * temp27 * temp3 * temp4 * temp33 * temp22;
              res28 += temp1 * temp2 * temp3 * temp4 * temp5 * temp34 * temp35;
              res29 += temp1 * temp8 * temp3 * temp4 * temp5 * temp36 * temp35;
              res30 += temp1 * temp2 * temp10 * temp4 * temp5 * temp37 * temp35;
              res31 += temp1 * temp8 * temp12 * temp4 * temp5 * temp37 * temp35;
              res32 += temp1 * temp13 * temp2 * temp14 * temp5 * temp37 * temp35;
              res33 += temp1 * temp2 * temp3 * temp15 * temp5 * temp38 * temp35;
              res34 += temp1 * temp8 * temp3 * temp17 * temp5 * temp38 * temp35;
              res35 += temp1 * temp2 * temp10 * temp18 * temp5 * temp38 * temp35;
              res36 += temp1 * temp8 * temp12 * temp18 * temp5 * temp38 * temp35;
              res37 += temp1 * temp8 * temp27 * temp3 * temp18 * temp30 * temp22;
              res38 += temp1 * temp2 * temp19 * temp3 * temp20 * temp38 * temp35;
              res39 += temp1 * temp8 * temp14 * temp3 * temp20 * temp38 * temp35;
              res40 += temp1 * temp2 * temp3 * temp4 * temp21 * temp39 * temp35;
              res41 += temp1 * temp8 * temp3 * temp4 * temp23 * temp39 * temp35;
              res42 += temp1 * temp2 * temp10 * temp4 * temp24 * temp39 * temp35;
              res43 += temp1 * temp8 * temp12 * temp4 * temp24 * temp39 * temp35;
              res44 += temp1 * temp13 * temp2 * temp14 * temp24 * temp39 * temp35;
              res45 += temp1 * temp2 * temp3 * temp15 * temp25 * temp39 * temp35;
              res46 += temp1 * temp8 * temp3 * temp17 * temp25 * temp39 * temp35;
              res47 += temp1 * temp2 * temp10 * temp18 * temp25 * temp39 * temp35;
              res48 += temp1 * temp8 * temp12 * temp18 * temp25 * temp39 * temp35;
            }
          }
        }
      }
    }
  }

  return 60 * (res11 + res10) + 45 * (res23 + res7 + res5) + 36 * (res13) + 24 * (res15 + res14) + 20 * (res27 + res4 + res3) + 18 * (res18 + res16) + 15 * (res34 + res32 + res25 + res24 + res22 + res21 + res9 + res8 + res6 + res2) + 10 * (res39 + res38 + res37 + res31 + res30 + res29 + res26) + 6 * (res41 + res20 + res19 + res17 + res12) + 5 * (res36 + res35 + res33) + 4 * (res43 + res42) + 3 * (res46 + res44) + res48 + res47 + res45 + res40 + res28 + res1;
}

case 8:
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  Affine2 res5(0.0);
  Affine2 res6(0.0);
  Affine2 res7(0.0);
  Affine2 res8(0.0);
  Affine2 res9(0.0);
  Affine2 res10(0.0);
  Affine2 res11(0.0);
  Affine2 res12(0.0);
  Affine2 res13(0.0);
  Affine2 res14(0.0);
  Affine2 res15(0.0);
  Affine2 res16(0.0);
  Affine2 res17(0.0);
  Affine2 res18(0.0);
  Affine2 res19(0.0);
  Affine2 res20(0.0);
  Affine2 res21(0.0);
  Affine2 res22(0.0);
  Affine2 res23(0.0);
  Affine2 res24(0.0);
  Affine2 res25(0.0);
  Affine2 res26(0.0);
  Affine2 res27(0.0);
  Affine2 res28(0.0);
  Affine2 res29(0.0);
  Affine2 res30(0.0);
  Affine2 res31(0.0);
  Affine2 res32(0.0);
  Affine2 res33(0.0);
  Affine2 res34(0.0);
  Affine2 res35(0.0);
  Affine2 res36(0.0);
  Affine2 res37(0.0);
  Affine2 res38(0.0);
  Affine2 res39(0.0);
  Affine2 res40(0.0);
  Affine2 res41(0.0);
  Affine2 res42(0.0);
  Affine2 res43(0.0);
  Affine2 res44(0.0);
  Affine2 res45(0.0);
  Affine2 res46(0.0);
  Affine2 res47(0.0);
  Affine2 res48(0.0);
  Affine2 res49(0.0);
  Affine2 res50(0.0);
  Affine2 res51(0.0);
  Affine2 res52(0.0);
  Affine2 res53(0.0);
  Affine2 res54(0.0);
  Affine2 res55(0.0);
  Affine2 res56(0.0);
  Affine2 res57(0.0);
  Affine2 res58(0.0);
  Affine2 res59(0.0);
  Affine2 res60(0.0);
  Affine2 res61(0.0);
  Affine2 res62(0.0);
  Affine2 res63(0.0);
  Affine2 res64(0.0);
  Affine2 res65(0.0);
  Affine2 res66(0.0);
  Affine2 res67(0.0);
  Affine2 res68(0.0);
  Affine2 res69(0.0);
  Affine2 res70(0.0);
  Affine2 res71(0.0);
  Affine2 res72(0.0);
  Affine2 res73(0.0);
  Affine2 res74(0.0);
  Affine2 res75(0.0);
  Affine2 res76(0.0);
  Affine2 res77(0.0);
  Affine2 res78(0.0);
  Affine2 res79(0.0);
  Affine2 res80(0.0);
  Affine2 res81(0.0);
  Affine2 res82(0.0);
  Affine2 res83(0.0);
  Affine2 res84(0.0);
  Affine2 res85(0.0);
  Affine2 res86(0.0);
  Affine2 res87(0.0);
  Affine2 res88(0.0);
  Affine2 res89(0.0);
  Affine2 res90(0.0);
  Affine2 res91(0.0);
  Affine2 res92(0.0);
  Affine2 res93(0.0);
  Affine2 res94(0.0);
  Affine2 res95(0.0);
  Affine2 res96(0.0);
  Affine2 res97(0.0);
  Affine2 res98(0.0);
  Affine2 res99(0.0);
  Affine2 res100(0.0);
  Affine2 res101(0.0);
  Affine2 res102(0.0);
  Affine2 res103(0.0);
  Affine2 res104(0.0);
  Affine2 res105(0.0);
  Affine2 res106(0.0);
  Affine2 res107(0.0);
  Affine2 res108(0.0);
  Affine2 res109(0.0);
  Affine2 res110(0.0);
  Affine2 res111(0.0);
  Affine2 res112(0.0);
  Affine2 res113(0.0);
  Affine2 res114(0.0);
  Affine2 res115(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
 std::vector<int> key3;
 std::vector<int> key4;
 std::vector<int> key5;
 std::vector<int> key6;
 std::vector<int> key7;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp7 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp51 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp6 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp36 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp56 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp5 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp23 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp55 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j3);
        Affine2 temp43 = edfr->eval_frechet(1, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp40 = edfr->eval_frechet(1, key0, y);

        for (int j4 = 0; j4 < nbvar; ++j4) {
          // Common Factor

          key0.clear();
          key0.push_back (j4);
          Affine2 temp4 = edfr->eval_frechet(0, key0, y);

          key0.clear();
          key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp17 = edfr->eval_frechet(4, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp54 = edfr->eval_frechet(3, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp46 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j4);
          Affine2 temp41 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp39 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j4);
          Affine2 temp31 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j3); key0.push_back (j4);
          Affine2 temp26 = edfr->eval_frechet(1, key0, y);

          for (int j5 = 0; j5 < nbvar; ++j5) {
            // Common Factor

            key0.clear();
            key0.push_back (j5);
            Affine2 temp3 = edfr->eval_frechet(0, key0, y);

            key0.clear();
            key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp12 = edfr->eval_frechet(5, key0, y);

            key0.clear();
            key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp53 = edfr->eval_frechet(4, key0, y);

            key0.clear();
            key0.push_back (j1); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp49 = edfr->eval_frechet(3, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp38 = edfr->eval_frechet(3, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp34 = edfr->eval_frechet(2, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j5);
            Affine2 temp29 = edfr->eval_frechet(1, key0, y);

            key0.clear();
            key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp25 = edfr->eval_frechet(2, key0, y);

            key0.clear();
            key0.push_back (j3); key0.push_back (j5);
            Affine2 temp21 = edfr->eval_frechet(1, key0, y);

            key0.clear();
            key0.push_back (j4); key0.push_back (j5);
            Affine2 temp19 = edfr->eval_frechet(1, key0, y);

            for (int j6 = 0; j6 < nbvar; ++j6) {
              // Common Factor

              key0.clear();
              key0.push_back (j6);
              Affine2 temp2 = edfr->eval_frechet(0, key0, y);

              key0.clear();
              key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp10 = edfr->eval_frechet(6, key0, y);

              key0.clear();
              key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp52 = edfr->eval_frechet(5, key0, y);

              key0.clear();
              key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp37 = edfr->eval_frechet(4, key0, y);

              key0.clear();
              key0.push_back (j2); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp44 = edfr->eval_frechet(3, key0, y);

              key0.clear();
              key0.push_back (j2); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp47 = edfr->eval_frechet(2, key0, y);

              key0.clear();
              key0.push_back (j2); key0.push_back (j6);
              Affine2 temp48 = edfr->eval_frechet(1, key0, y);

              key0.clear();
              key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp24 = edfr->eval_frechet(3, key0, y);

              key0.clear();
              key0.push_back (j3); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp32 = edfr->eval_frechet(2, key0, y);

              key0.clear();
              key0.push_back (j3); key0.push_back (j6);
              Affine2 temp28 = edfr->eval_frechet(1, key0, y);

              key0.clear();
              key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp18 = edfr->eval_frechet(2, key0, y);

              key0.clear();
              key0.push_back (j4); key0.push_back (j6);
              Affine2 temp15 = edfr->eval_frechet(1, key0, y);

              key0.clear();
              key0.push_back (j5); key0.push_back (j6);
              Affine2 temp13 = edfr->eval_frechet(1, key0, y);

              for (int j7 = 0; j7 < nbvar; ++j7) {
                // Common Factor

                key0.clear();
                key0.push_back (j7);
                Affine2 temp1 = edfr->eval_frechet(0, key0, y);

                key0.clear();
                key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp8 = edfr->eval_frechet(7, key0, y);

                key0.clear();
                key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp50 = edfr->eval_frechet(6, key0, y);

                key0.clear();
                key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp35 = edfr->eval_frechet(5, key0, y);

                key0.clear();
                key0.push_back (j2); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp42 = edfr->eval_frechet(4, key0, y);

                key0.clear();
                key0.push_back (j2); key0.push_back (j5); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp45 = edfr->eval_frechet(3, key0, y);

                key0.clear();
                key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp22 = edfr->eval_frechet(4, key0, y);

                key0.clear();
                key0.push_back (j3); key0.push_back (j5); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp30 = edfr->eval_frechet(3, key0, y);

                key0.clear();
                key0.push_back (j3); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp33 = edfr->eval_frechet(2, key0, y);

                key0.clear();
                key0.push_back (j4); key0.push_back (j5); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp16 = edfr->eval_frechet(3, key0, y);

                key0.clear();
                key0.push_back (j4); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp20 = edfr->eval_frechet(2, key0, y);

                key0.clear();
                key0.push_back (j4); key0.push_back (j7);
                Affine2 temp27 = edfr->eval_frechet(1, key0, y);

                key0.clear();
                key0.push_back (j5); key0.push_back (j6); key0.push_back (j7);
                Affine2 temp11 = edfr->eval_frechet(2, key0, y);

                key0.clear();
                key0.push_back (j5); key0.push_back (j7);
                Affine2 temp14 = edfr->eval_frechet(1, key0, y);

                key0.clear();
                key0.push_back (j6); key0.push_back (j7);
                Affine2 temp9 = edfr->eval_frechet(1, key0, y);

                // Tree computation
                res1 += temp1 * temp2 * temp3 * temp4 * temp5 * temp6 * temp7 * temp8;
                res2 += temp1 * temp9 * temp3 * temp4 * temp5 * temp6 * temp7 * temp10;
                res3 += temp1 * temp2 * temp11 * temp4 * temp5 * temp6 * temp7 * temp12;
                res4 += temp1 * temp9 * temp13 * temp4 * temp5 * temp6 * temp7 * temp12;
                res5 += temp1 * temp14 * temp2 * temp15 * temp5 * temp6 * temp7 * temp12;
                res6 += temp1 * temp2 * temp3 * temp16 * temp5 * temp6 * temp7 * temp17;
                res7 += temp1 * temp9 * temp3 * temp18 * temp5 * temp6 * temp7 * temp17;
                res8 += temp1 * temp2 * temp11 * temp19 * temp5 * temp6 * temp7 * temp17;
                res9 += temp1 * temp9 * temp13 * temp19 * temp5 * temp6 * temp7 * temp17;
                res10 += temp1 * temp2 * temp20 * temp3 * temp21 * temp6 * temp7 * temp17;
                res11 += temp1 * temp9 * temp15 * temp3 * temp21 * temp6 * temp7 * temp17;
                res12 += temp1 * temp2 * temp3 * temp4 * temp22 * temp6 * temp7 * temp23;
                res13 += temp1 * temp9 * temp3 * temp4 * temp24 * temp6 * temp7 * temp23;
                res14 += temp1 * temp2 * temp11 * temp4 * temp25 * temp6 * temp7 * temp23;
                res15 += temp1 * temp9 * temp13 * temp4 * temp25 * temp6 * temp7 * temp23;
                res16 += temp1 * temp14 * temp2 * temp15 * temp25 * temp6 * temp7 * temp23;
                res17 += temp1 * temp2 * temp3 * temp16 * temp26 * temp6 * temp7 * temp23;
                res18 += temp1 * temp9 * temp3 * temp18 * temp26 * temp6 * temp7 * temp23;
                res19 += temp1 * temp2 * temp11 * temp19 * temp26 * temp6 * temp7 * temp23;
                res20 += temp1 * temp9 * temp13 * temp19 * temp26 * temp6 * temp7 * temp23;
                res21 += temp1 * temp27 * temp2 * temp28 * temp3 * temp29 * temp7 * temp17;
                res22 += temp1 * temp2 * temp3 * temp30 * temp4 * temp31 * temp7 * temp23;
                res23 += temp1 * temp9 * temp3 * temp32 * temp4 * temp31 * temp7 * temp23;
                res24 += temp1 * temp2 * temp11 * temp21 * temp4 * temp31 * temp7 * temp23;
                res25 += temp1 * temp9 * temp13 * temp21 * temp4 * temp31 * temp7 * temp23;
                res26 += temp1 * temp2 * temp33 * temp3 * temp4 * temp34 * temp7 * temp23;
                res27 += temp1 * temp9 * temp28 * temp3 * temp4 * temp34 * temp7 * temp23;
                res28 += temp1 * temp2 * temp3 * temp4 * temp5 * temp35 * temp7 * temp36;
                res29 += temp1 * temp9 * temp3 * temp4 * temp5 * temp37 * temp7 * temp36;
                res30 += temp1 * temp2 * temp11 * temp4 * temp5 * temp38 * temp7 * temp36;
                res31 += temp1 * temp9 * temp13 * temp4 * temp5 * temp38 * temp7 * temp36;
                res32 += temp1 * temp14 * temp2 * temp15 * temp5 * temp38 * temp7 * temp36;
                res33 += temp1 * temp2 * temp3 * temp16 * temp5 * temp39 * temp7 * temp36;
                res34 += temp1 * temp9 * temp3 * temp18 * temp5 * temp39 * temp7 * temp36;
                res35 += temp1 * temp2 * temp11 * temp19 * temp5 * temp39 * temp7 * temp36;
                res36 += temp1 * temp9 * temp13 * temp19 * temp5 * temp39 * temp7 * temp36;
                res37 += temp1 * temp9 * temp28 * temp3 * temp19 * temp31 * temp7 * temp23;
                res38 += temp1 * temp2 * temp20 * temp3 * temp21 * temp39 * temp7 * temp36;
                res39 += temp1 * temp9 * temp15 * temp3 * temp21 * temp39 * temp7 * temp36;
                res40 += temp1 * temp2 * temp3 * temp4 * temp22 * temp40 * temp7 * temp36;
                res41 += temp1 * temp9 * temp3 * temp4 * temp24 * temp40 * temp7 * temp36;
                res42 += temp1 * temp2 * temp11 * temp4 * temp25 * temp40 * temp7 * temp36;
                res43 += temp1 * temp9 * temp13 * temp4 * temp25 * temp40 * temp7 * temp36;
                res44 += temp1 * temp14 * temp2 * temp15 * temp25 * temp40 * temp7 * temp36;
                res45 += temp1 * temp2 * temp3 * temp16 * temp26 * temp40 * temp7 * temp36;
                res46 += temp1 * temp9 * temp3 * temp18 * temp26 * temp40 * temp7 * temp36;
                res47 += temp1 * temp2 * temp11 * temp19 * temp26 * temp40 * temp7 * temp36;
                res48 += temp1 * temp9 * temp13 * temp19 * temp26 * temp40 * temp7 * temp36;
                res49 += temp1 * temp2 * temp33 * temp3 * temp29 * temp4 * temp41 * temp23;
                res50 += temp1 * temp9 * temp28 * temp3 * temp29 * temp4 * temp41 * temp23;
                res51 += temp1 * temp2 * temp3 * temp4 * temp42 * temp5 * temp43 * temp36;
                res52 += temp1 * temp9 * temp3 * temp4 * temp44 * temp5 * temp43 * temp36;
                res53 += temp1 * temp2 * temp11 * temp4 * temp34 * temp5 * temp43 * temp36;
                res54 += temp1 * temp9 * temp13 * temp4 * temp34 * temp5 * temp43 * temp36;
                res55 += temp1 * temp14 * temp2 * temp15 * temp34 * temp5 * temp43 * temp36;
                res56 += temp1 * temp2 * temp3 * temp16 * temp31 * temp5 * temp43 * temp36;
                res57 += temp1 * temp9 * temp3 * temp18 * temp31 * temp5 * temp43 * temp36;
                res58 += temp1 * temp2 * temp11 * temp19 * temp31 * temp5 * temp43 * temp36;
                res59 += temp1 * temp9 * temp13 * temp19 * temp31 * temp5 * temp43 * temp36;
                res60 += temp1 * temp2 * temp3 * temp45 * temp4 * temp5 * temp46 * temp36;
                res61 += temp1 * temp9 * temp3 * temp47 * temp4 * temp5 * temp46 * temp36;
                res62 += temp1 * temp2 * temp11 * temp29 * temp4 * temp5 * temp46 * temp36;
                res63 += temp1 * temp9 * temp13 * temp29 * temp4 * temp5 * temp46 * temp36;
                res64 += temp1 * temp9 * temp48 * temp3 * temp4 * temp5 * temp49 * temp36;
                res65 += temp1 * temp2 * temp3 * temp4 * temp5 * temp6 * temp50 * temp51;
                res66 += temp1 * temp9 * temp3 * temp4 * temp5 * temp6 * temp52 * temp51;
                res67 += temp1 * temp2 * temp11 * temp4 * temp5 * temp6 * temp53 * temp51;
                res68 += temp1 * temp9 * temp13 * temp4 * temp5 * temp6 * temp53 * temp51;
                res69 += temp1 * temp14 * temp2 * temp15 * temp5 * temp6 * temp53 * temp51;
                res70 += temp1 * temp2 * temp3 * temp16 * temp5 * temp6 * temp54 * temp51;
                res71 += temp1 * temp9 * temp3 * temp18 * temp5 * temp6 * temp54 * temp51;
                res72 += temp1 * temp2 * temp11 * temp19 * temp5 * temp6 * temp54 * temp51;
                res73 += temp1 * temp9 * temp13 * temp19 * temp5 * temp6 * temp54 * temp51;
                res74 += temp1 * temp9 * temp48 * temp3 * temp19 * temp5 * temp46 * temp36;
                res75 += temp1 * temp2 * temp20 * temp3 * temp21 * temp6 * temp54 * temp51;
                res76 += temp1 * temp9 * temp15 * temp3 * temp21 * temp6 * temp54 * temp51;
                res77 += temp1 * temp2 * temp3 * temp4 * temp22 * temp6 * temp55 * temp51;
                res78 += temp1 * temp9 * temp3 * temp4 * temp24 * temp6 * temp55 * temp51;
                res79 += temp1 * temp2 * temp11 * temp4 * temp25 * temp6 * temp55 * temp51;
                res80 += temp1 * temp9 * temp13 * temp4 * temp25 * temp6 * temp55 * temp51;
                res81 += temp1 * temp14 * temp2 * temp15 * temp25 * temp6 * temp55 * temp51;
                res82 += temp1 * temp2 * temp3 * temp16 * temp26 * temp6 * temp55 * temp51;
                res83 += temp1 * temp9 * temp3 * temp18 * temp26 * temp6 * temp55 * temp51;
                res84 += temp1 * temp2 * temp11 * temp19 * temp26 * temp6 * temp55 * temp51;
                res85 += temp1 * temp9 * temp13 * temp19 * temp26 * temp6 * temp55 * temp51;
                res86 += temp1 * temp2 * temp11 * temp29 * temp4 * temp26 * temp43 * temp36;
                res87 += temp1 * temp9 * temp13 * temp29 * temp4 * temp26 * temp43 * temp36;
                res88 += temp1 * temp27 * temp2 * temp28 * temp3 * temp29 * temp54 * temp51;
                res89 += temp1 * temp2 * temp3 * temp30 * temp4 * temp31 * temp55 * temp51;
                res90 += temp1 * temp9 * temp3 * temp32 * temp4 * temp31 * temp55 * temp51;
                res91 += temp1 * temp2 * temp11 * temp21 * temp4 * temp31 * temp55 * temp51;
                res92 += temp1 * temp9 * temp13 * temp21 * temp4 * temp31 * temp55 * temp51;
                res93 += temp1 * temp2 * temp33 * temp3 * temp4 * temp34 * temp55 * temp51;
                res94 += temp1 * temp9 * temp28 * temp3 * temp4 * temp34 * temp55 * temp51;
                res95 += temp1 * temp2 * temp3 * temp4 * temp5 * temp35 * temp56 * temp51;
                res96 += temp1 * temp9 * temp3 * temp4 * temp5 * temp37 * temp56 * temp51;
                res97 += temp1 * temp2 * temp11 * temp4 * temp5 * temp38 * temp56 * temp51;
                res98 += temp1 * temp9 * temp13 * temp4 * temp5 * temp38 * temp56 * temp51;
                res99 += temp1 * temp14 * temp2 * temp15 * temp5 * temp38 * temp56 * temp51;
                res100 += temp1 * temp2 * temp3 * temp16 * temp5 * temp39 * temp56 * temp51;
                res101 += temp1 * temp9 * temp3 * temp18 * temp5 * temp39 * temp56 * temp51;
                res102 += temp1 * temp2 * temp11 * temp19 * temp5 * temp39 * temp56 * temp51;
                res103 += temp1 * temp9 * temp13 * temp19 * temp5 * temp39 * temp56 * temp51;
                res104 += temp1 * temp9 * temp28 * temp3 * temp19 * temp31 * temp55 * temp51;
                res105 += temp1 * temp2 * temp20 * temp3 * temp21 * temp39 * temp56 * temp51;
                res106 += temp1 * temp9 * temp15 * temp3 * temp21 * temp39 * temp56 * temp51;
                res107 += temp1 * temp2 * temp3 * temp4 * temp22 * temp40 * temp56 * temp51;
                res108 += temp1 * temp9 * temp3 * temp4 * temp24 * temp40 * temp56 * temp51;
                res109 += temp1 * temp2 * temp11 * temp4 * temp25 * temp40 * temp56 * temp51;
                res110 += temp1 * temp9 * temp13 * temp4 * temp25 * temp40 * temp56 * temp51;
                res111 += temp1 * temp14 * temp2 * temp15 * temp25 * temp40 * temp56 * temp51;
                res112 += temp1 * temp2 * temp3 * temp16 * temp26 * temp40 * temp56 * temp51;
                res113 += temp1 * temp9 * temp3 * temp18 * temp26 * temp40 * temp56 * temp51;
                res114 += temp1 * temp2 * temp11 * temp19 * temp26 * temp40 * temp56 * temp51;
                res115 += temp1 * temp9 * temp13 * temp19 * temp26 * temp40 * temp56 * temp51;
              }
            }
          }
        }
      }
    }
  }

  return 315 * (res23) + 210 * (res11 + res10) + 140 * (res27) + 126 * (res52 + res13) + 105 * (res74 + res61 + res50 + res49 + res34 + res32 + res25 + res24 + res22 + res21 + res7 + res5) + 84 * (res54 + res53 + res15 + res14) + 70 * (res39 + res38 + res37 + res31 + res30 + res29 + res26) + 63 * (res57 + res55 + res18 + res16) + 60 * (res76 + res75) + 45 * (res90 + res71 + res69) + 42 * (res41) + 36 * (res78) + 35 * (res87 + res86 + res64 + res63 + res62 + res60 + res36 + res35 + res33 + res9 + res8 + res6 + res4 + res3) + 28 * (res43 + res42) + 24 * (res80 + res79) + 21 * (res59 + res58 + res56 + res51 + res46 + res44 + res20 + res19 + res17 + res12 + res2) + 20 * (res94 + res68 + res67) + 18 * (res83 + res81) + 15 * (res101 + res99 + res92 + res91 + res89 + res88 + res73 + res72 + res70 + res66) + 10 * (res106 + res105 + res104 + res98 + res97 + res96 + res93) + 7 * (res48 + res47 + res45 + res40 + res28) + 6 * (res108 + res85 + res84 + res82 + res77) + 5 * (res103 + res102 + res100) + 4 * (res110 + res109) + 3 * (res113 + res111) + res115 + res114 + res112 + res107 + res95 + res65 + res1;
}



default:
{
//std::cerr<< "Not implemented" << std::endl;
return Affine2(Interval::ALL_REALS);
}
    }
};


// Butcher table of ExplicitRK4 with the form:
// c | A
// ------
//   | b
// is defined by
// c = [ 0; 1/2; 1/2; 1 ]
// A = [ [ 0; 0; 0; 0 ];   [ 1/2; 0; 0; 0 ];   [ 0; 1/2; 0; 0 ];   [ 0; 0; 1; 0 ] ]
// b = [ 1/6; 1/3; 1/3; 1/6 ]


Affine2 edtree_frechet::lteExplicitRK4 (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  Affine2 res5(0.0);
  Affine2 res6(0.0);
  Affine2 res7(0.0);
  Affine2 res8(0.0);
  Affine2 res9(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
 std::vector<int> key3;
 std::vector<int> key4;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp4 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp14 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp3 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp9 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp16 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp2 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp7 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp15 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j3);
        Affine2 temp12 = edfr->eval_frechet(1, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp10 = edfr->eval_frechet(1, key0, y);

        for (int j4 = 0; j4 < nbvar; ++j4) {
          // Common Factor

          key0.clear();
          key0.push_back (j4);
          Affine2 temp1 = edfr->eval_frechet(0, key0, y);

          key0.clear();
          key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp5 = edfr->eval_frechet(4, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp13 = edfr->eval_frechet(3, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp8 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j4);
          Affine2 temp11 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j3); key0.push_back (j4);
          Affine2 temp6 = edfr->eval_frechet(1, key0, y);

          // Tree computation
          res1 += temp1 * temp2 * temp3 * temp4 * temp5;
          res2 += temp1 * temp6 * temp3 * temp4 * temp7;
          res3 += temp1 * temp2 * temp8 * temp4 * temp9;
          res4 += temp1 * temp6 * temp10 * temp4 * temp9;
          res5 += temp1 * temp11 * temp2 * temp12 * temp9;
          res6 += temp1 * temp2 * temp3 * temp13 * temp14;
          res7 += temp1 * temp6 * temp3 * temp15 * temp14;
          res8 += temp1 * temp2 * temp8 * temp16 * temp14;
          res9 += temp1 * temp6 * temp10 * temp16 * temp14;
        }
      }
    }
  }

  return res9 + ((double)1/2) * (res7) + ((double)1/4) * (res3) + ((double)1/6) * (res6) + ((double)-1/24) * (res1) + ((double)-1/4) * (res8 + res2) + ((double)-3/4) * (res5) + ((double)-1) * (res4);
}


// Butcher table of ImplicitLobbato3a4 with the form:
// c | A
// ------
//   | b
// is defined by
// c = [ 0; 1/2; 1 ]
// A = [ [ 0; 0; 0 ];   [ 5/24; 1/3; -1/24 ];   [ 1/6; 2/3; 1/6 ] ]
// b = [ 1/6; 2/3; 1/6 ]


Affine2 edtree_frechet::lteImplicitLobbato3a4 (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  Affine2 res5(0.0);
  Affine2 res6(0.0);
  Affine2 res7(0.0);
  Affine2 res8(0.0);
  Affine2 res9(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
 std::vector<int> key3;
 std::vector<int> key4;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp4 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp14 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp3 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp9 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp16 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp2 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp7 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp15 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j3);
        Affine2 temp12 = edfr->eval_frechet(1, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp10 = edfr->eval_frechet(1, key0, y);

        for (int j4 = 0; j4 < nbvar; ++j4) {
          // Common Factor

          key0.clear();
          key0.push_back (j4);
          Affine2 temp1 = edfr->eval_frechet(0, key0, y);

          key0.clear();
          key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp5 = edfr->eval_frechet(4, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp13 = edfr->eval_frechet(3, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp8 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j4);
          Affine2 temp11 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j3); key0.push_back (j4);
          Affine2 temp6 = edfr->eval_frechet(1, key0, y);

          // Tree computation
          res1 += temp1 * temp2 * temp3 * temp4 * temp5;
          res2 += temp1 * temp6 * temp3 * temp4 * temp7;
          res3 += temp1 * temp2 * temp8 * temp4 * temp9;
          res4 += temp1 * temp6 * temp10 * temp4 * temp9;
          res5 += temp1 * temp11 * temp2 * temp12 * temp9;
          res6 += temp1 * temp2 * temp3 * temp13 * temp14;
          res7 += temp1 * temp6 * temp3 * temp15 * temp14;
          res8 += temp1 * temp2 * temp8 * temp16 * temp14;
          res9 += temp1 * temp6 * temp10 * temp16 * temp14;
        }
      }
    }
  }

  return ((double)1.0/2.0) * (res7) + ((double)1.0/6.0) * (res9 + res8 + res6) + ((double)-1.0/24.0) * (res1) + ((double)-1.0/8.0) * (res5) + ((double)-1.0/6.0) * (res4 + res3) + ((double)-1.0/4.0) * (res2);
}



// Butcher table of ImplicitEuler with the form:
// c | A
// ------
//   | b
// is defined by
// c = [ 1 ]
// A = [ [ 1 ] ]
// b = [ 1 ]


Affine2 edtree_frechet::lteImplicitEuler (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp1 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp2 = edfr->eval_frechet(1, key0, y);

    // Tree computation
    res1 += temp1 * temp2;
  }

  return ((double)-1) * (res1);
}

// Butcher table of ImplicitMidpoint with the form:
// c | A
// ------
//   | b
// is defined by
// c = [ 1/2 ]
// A = [ [ 1/2 ] ]
// b = [ 1 ]


Affine2 edtree_frechet::lteImplicitMidpoint (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp2 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp5 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp1 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp3 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp4 = edfr->eval_frechet(1, key0, y);

      // Tree computation
      res1 += temp1 * temp2 * temp3;
      res2 += temp1 * temp4 * temp5;
    }
  }

  return ((double)1/4) * (res1) + ((double)-1/2) * (res2);
}

// Butcher table of ImplicitRadau3 with the form:
// c | A
// ------
//   | b
// is defined by
// c = [ 1/3; 1 ]
// A = [ [ 5/12; -1/12 ];   [ 3/4; 1/4 ] ]
// b = [ 3/4; 1/4 ]


Affine2 edtree_frechet::lteImplicitRadau3 (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
 std::vector<int> key3;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp3 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp8 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp2 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp6 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp9 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp1 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp4 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp7 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp5 = edfr->eval_frechet(1, key0, y);

        // Tree computation
        res1 += temp1 * temp2 * temp3 * temp4;
        res2 += temp1 * temp5 * temp3 * temp6;
        res3 += temp1 * temp2 * temp7 * temp8;
        res4 += temp1 * temp5 * temp9 * temp8;
      }
    }
  }

  return (1.0/3.0) * (res4 + res3) + (-1.0/9.0) * (res1) + (-1.0/3.0) * (res2);
}



// Butcher table of ImplicitLobbato3c4 with the form:
// c | A
// ------
//   | b
// is defined by
// c = [ 0; 1/2; 1 ]
// A = [ [ 1/6; -1/3; 1/6 ];   [ 1/6; 5/12; -1/12 ];   [ 1/6; 2/3; 1/6 ] ]
// b = [ 1/6; 2/3; 1/6 ]


Affine2 edtree_frechet::lteImplicitLobbato3c4 (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  Affine2 res5(0.0);
  Affine2 res6(0.0);
  Affine2 res7(0.0);
  Affine2 res8(0.0);
  Affine2 res9(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
 std::vector<int> key3;
 std::vector<int> key4;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp4 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp14 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp3 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp9 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp16 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp2 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp7 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp15 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j3);
        Affine2 temp12 = edfr->eval_frechet(1, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp10 = edfr->eval_frechet(1, key0, y);

        for (int j4 = 0; j4 < nbvar; ++j4) {
          // Common Factor

          key0.clear();
          key0.push_back (j4);
          Affine2 temp1 = edfr->eval_frechet(0, key0, y);

          key0.clear();
          key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp5 = edfr->eval_frechet(4, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp13 = edfr->eval_frechet(3, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp8 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j4);
          Affine2 temp11 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j3); key0.push_back (j4);
          Affine2 temp6 = edfr->eval_frechet(1, key0, y);

          // Tree computation
          res1 += temp1 * temp2 * temp3 * temp4 * temp5;
          res2 += temp1 * temp6 * temp3 * temp4 * temp7;
          res3 += temp1 * temp2 * temp8 * temp4 * temp9;
          res4 += temp1 * temp6 * temp10 * temp4 * temp9;
          res5 += temp1 * temp11 * temp2 * temp12 * temp9;
          res6 += temp1 * temp2 * temp3 * temp13 * temp14;
          res7 += temp1 * temp6 * temp3 * temp15 * temp14;
          res8 += temp1 * temp2 * temp8 * temp16 * temp14;
          res9 += temp1 * temp6 * temp10 * temp16 * temp14;
        }
      }
    }
  }

  return ((double)1./2) * (res7) + ((double)1./4) * (res4 + res3) + ((double)1./6) * (res6) + ((double)-1./24) * (res1) + ((double)-1./8) * (res5) + ((double)-1./4) * (res9 + res8 + res2);
}

// Butcher table of ExplicitHeun with the form:
// c | A
// ------
//   | b
// is defined by
// c = [ 0; 1 ]
// A = [ [ 0; 0 ];   [ 1; 0 ] ]
// b = [ 1/2; 1/2 ]


Affine2 edtree_frechet::lteExplicitHeun (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp2 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp5 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp1 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp3 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp4 = edfr->eval_frechet(1, key0, y);

      // Tree computation
      res1 += temp1 * temp2 * temp3;
      res2 += temp1 * temp4 * temp5;
    }
  }

  return res2 + ((double)-1.0/2.0) * (res1);
}

Affine2 edtree_frechet::lteImplicitRadau2a5 (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  Affine2 res5(0.0);
  Affine2 res6(0.0);
  Affine2 res7(0.0);
  Affine2 res8(0.0);
  Affine2 res9(0.0);
  Affine2 res10(0.0);
  Affine2 res11(0.0);
  Affine2 res12(0.0);
  Affine2 res13(0.0);
  Affine2 res14(0.0);
  Affine2 res15(0.0);
  Affine2 res16(0.0);
  Affine2 res17(0.0);
  Affine2 res18(0.0);
  Affine2 res19(0.0);
  Affine2 res20(0.0);
 std::vector<int> key0;
 std::vector<int> key1;
 std::vector<int> key2;
 std::vector<int> key3;
 std::vector<int> key4;
 std::vector<int> key5;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp5 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp21 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp4 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp15 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp24 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp3 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp10 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp23 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j3);
        Affine2 temp19 = edfr->eval_frechet(1, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp17 = edfr->eval_frechet(1, key0, y);

        for (int j4 = 0; j4 < nbvar; ++j4) {
          // Common Factor

          key0.clear();
          key0.push_back (j4);
          Affine2 temp2 = edfr->eval_frechet(0, key0, y);

          key0.clear();
          key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp8 = edfr->eval_frechet(4, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp22 = edfr->eval_frechet(3, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp16 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j4);
          Affine2 temp13 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j3); key0.push_back (j4);
          Affine2 temp11 = edfr->eval_frechet(1, key0, y);

          for (int j5 = 0; j5 < nbvar; ++j5) {
            // Common Factor

            key0.clear();
            key0.push_back (j5);
            Affine2 temp1 = edfr->eval_frechet(0, key0, y);

            key0.clear();
            key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp6 = edfr->eval_frechet(5, key0, y);

            key0.clear();
            key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp20 = edfr->eval_frechet(4, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp14 = edfr->eval_frechet(3, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp18 = edfr->eval_frechet(2, key0, y);

            key0.clear();
            key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp9 = edfr->eval_frechet(2, key0, y);

            key0.clear();
            key0.push_back (j3); key0.push_back (j5);
            Affine2 temp12 = edfr->eval_frechet(1, key0, y);

            key0.clear();
            key0.push_back (j4); key0.push_back (j5);
            Affine2 temp7 = edfr->eval_frechet(1, key0, y);

            // Tree computation
            res1 += temp1 * temp2 * temp3 * temp4 * temp5 * temp6;
            res2 += temp1 * temp7 * temp3 * temp4 * temp5 * temp8;
            res3 += temp1 * temp2 * temp9 * temp4 * temp5 * temp10;
            res4 += temp1 * temp7 * temp11 * temp4 * temp5 * temp10;
            res5 += temp1 * temp12 * temp2 * temp13 * temp5 * temp10;
            res6 += temp1 * temp2 * temp3 * temp14 * temp5 * temp15;
            res7 += temp1 * temp7 * temp3 * temp16 * temp5 * temp15;
            res8 += temp1 * temp2 * temp9 * temp17 * temp5 * temp15;
            res9 += temp1 * temp7 * temp11 * temp17 * temp5 * temp15;
            res10 += temp1 * temp2 * temp18 * temp3 * temp19 * temp15;
            res11 += temp1 * temp7 * temp13 * temp3 * temp19 * temp15;
            res12 += temp1 * temp2 * temp3 * temp4 * temp20 * temp21;
            res13 += temp1 * temp7 * temp3 * temp4 * temp22 * temp21;
            res14 += temp1 * temp2 * temp9 * temp4 * temp23 * temp21;
            res15 += temp1 * temp7 * temp11 * temp4 * temp23 * temp21;
            res16 += temp1 * temp12 * temp2 * temp13 * temp23 * temp21;
            res17 += temp1 * temp2 * temp3 * temp14 * temp24 * temp21;
            res18 += temp1 * temp7 * temp3 * temp16 * temp24 * temp21;
            res19 += temp1 * temp2 * temp9 * temp17 * temp24 * temp21;
            res20 += temp1 * temp7 * temp11 * temp17 * temp24 * temp21;
          }
        }
      }
    }
  }

  return Interval(-0.10000000000002673E-1, -0.99999999999984546E-2) * res1 +
  Interval(-0.10000000000002673E0, -0.99999999999982325E-1) * res2 +
  Interval(-0.10000000000002895E0, -0.99999999999980105E-1) * res3 +
  Interval(-0.10000000000003340E0, -0.99999999999973443E-1) * res4 +
  Interval(-0.15000000000004010E0, -0.14999999999997016E0) * res5 +
  Interval(0.99999999999983435E-1, 0.10000000000001286E0) * res6 +
  Interval(0.29999999999994198E0, 0.30000000000004357E0) * res7 +
  Interval(0.99999999999975664E-1, 0.10000000000002229E0) * res8 +
  Interval(0.99999999999966227E-1, 0.10000000000003229E0) * res9 +
  Interval(-0.10000000000002895E0, -0.99999999999980105E-1) * res10 +
  Interval(-0.10000000000003340E0, -0.99999999999971223E-1) * res11 +
  Interval(0.49999999999994826E-1, 0.50000000000004374E-1) * res12 +
  Interval(0.29999999999996829E0, 0.30000000000002758E0) * res13 +
  Interval(0.19999999999997531E0, 0.20000000000002238E0) * res14 +
  Interval(0.19999999999996643E0, 0.20000000000003082E0) * res15 +
  Interval(0.14999999999998281E0, 0.15000000000001479E0) * res16 +
  Interval(-0.10000000000001075E0, -0.99999999999990097E-1) * res17 +
  Interval(-0.30000000000003690E0, -0.29999999999996563E0) * res18 +
  Interval(-0.10000000000001741E0, -0.99999999999982547E-1) * res19 +
  Interval(-0.10000000000002451E0, -0.99999999999975220E-1) * res20;
}


// Butcher table of implicitGauss4 with the form:
// c | A
// ------
//   | b
// is defined by
// c = [ 8453/40000; 31547/40000 ]
// A = [ [ 1/4; -1547/40000 ];   [ 21547/40000; 1/4 ] ]
// b = [ 1/2; 1/2 ]


Affine2 edtree_frechet::lteImplicitGauss4 (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  Affine2 res5(0.0);
  Affine2 res6(0.0);
  Affine2 res7(0.0);
  Affine2 res8(0.0);
  Affine2 res9(0.0);
  std::vector<int> key0;
  std::vector<int> key1;
  std::vector<int> key2;
  std::vector<int> key3;
  std::vector<int> key4;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp4 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp14 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp3 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp9 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp16 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp2 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp7 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp15 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j3);
        Affine2 temp12 = edfr->eval_frechet(1, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp10 = edfr->eval_frechet(1, key0, y);

        for (int j4 = 0; j4 < nbvar; ++j4) {
          // Common Factor

          key0.clear();
          key0.push_back (j4);
          Affine2 temp1 = edfr->eval_frechet(0, key0, y);

          key0.clear();
          key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp5 = edfr->eval_frechet(4, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp13 = edfr->eval_frechet(3, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp8 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j4);
          Affine2 temp11 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j3); key0.push_back (j4);
          Affine2 temp6 = edfr->eval_frechet(1, key0, y);

          // Tree computation
          res1 += temp1 * temp2 * temp3 * temp4 * temp5;
          res2 += temp1 * temp6 * temp3 * temp4 * temp7;
          res3 += temp1 * temp2 * temp8 * temp4 * temp9;
          res4 += temp1 * temp6 * temp10 * temp4 * temp9;
          res5 += temp1 * temp11 * temp2 * temp12 * temp9;
          res6 += temp1 * temp2 * temp3 * temp13 * temp14;
          res7 += temp1 * temp6 * temp3 * temp15 * temp14;
          res8 += temp1 * temp2 * temp8 * temp16 * temp14;
          res9 += temp1 * temp6 * temp10 * temp16 * temp14;
        }
      }
    }
  }

  return ((double)21333532266713043/128000000000000000) * (res2) + ((double)21332935466713043/128000000000000000) * (res8) + ((double)10666318533286957/64000000000000000) * (res9) + ((double)10666616933286957/128000000000000000) * (res5) + ((double)14222553777762319/512000000000000000) * (res1) + ((double)-14222255377762319/128000000000000000) * (res6) + ((double)-21333233866713043/128000000000000000) * (res3) + ((double)-10666915333286957/64000000000000000) * (res4) + ((double)-21333532266713043/64000000000000000) * (res7);
}

Affine2 edtree_frechet::lteexplicitkutta3 (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  std::vector<int> key0;
  std::vector<int> key1;
  std::vector<int> key2;
  std::vector<int> key3;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp3 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp8 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

//       key0.clear();
//       key0.push_back (j2);
//       Affine2 temp2 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp6 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp9 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp1 = edfr->eval_frechet(0, key0, y);

//         key0.clear();
//         key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
//         Affine2 temp4 = edfr->eval_frechet(3, key0, y);
//
//         key0.clear();
//         key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
//         Affine2 temp7 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp5 = edfr->eval_frechet(1, key0, y);

        // Tree computation
        //res1 += temp1 * temp2 * temp3 * temp4;
        res2 += temp1 * temp5 * temp3 * temp6;
        //res3 += temp1 * temp2 * temp7 * temp8;
        res4 += temp1 * temp5 * temp9 * temp8;
      }
    }
  }

  return /*Interval(1) * res1 +*/ Interval(-1.0) * res2 + /*Interval(1) * res3 + Interval(1) **/ res4;
}

// Butcher table of implicitGauss6 with the form:
// c | A
// ------
//   | b
// is defined by
// c = [ Interval(0.11270166537925830E0, 0.11270166537925835E0); Interval(0.50000000000000000E0, 0.50000000000000000E0); Interval(0.88729833462074159E0, 0.88729833462074170E0) ]
// A = [ [ Interval(0.13888888888888887E0, 0.13888888888888890E0); Interval(-0.35976667524938943E-1, -0.35976667524938860E-1); Interval(0.97894440153082907E-2, 0.97894440153083462E-2) ];   [ Interval(0.30026319498086451E0, 0.30026319498086462E0); Interval(0.22222222222222221E0, 0.22222222222222224E0); Interval(-0.22485417203086860E-1, -0.22485417203086777E-1) ];   [ Interval(0.30026319498086451E0, 0.30026319498086462E0); Interval(0.48042111196938331E0, 0.48042111196938342E0); Interval(0.13888888888888887E0, 0.13888888888888890E0) ] ]
// b = [ Interval(0.27777777777777773E0, 0.27777777777777779E0); Interval(0.44444444444444442E0, 0.44444444444444448E0); Interval(0.27777777777777773E0, 0.27777777777777779E0) ]


Affine2 edtree_frechet::lteImplicitGauss6 (int j, Affine2Vector y)
{
  Affine2 res1(0.0);
  Affine2 res2(0.0);
  Affine2 res3(0.0);
  Affine2 res4(0.0);
  Affine2 res5(0.0);
  Affine2 res6(0.0);
  Affine2 res7(0.0);
  Affine2 res8(0.0);
  Affine2 res9(0.0);
  Affine2 res10(0.0);
  Affine2 res11(0.0);
  Affine2 res12(0.0);
  Affine2 res13(0.0);
  Affine2 res14(0.0);
  Affine2 res15(0.0);
  Affine2 res16(0.0);
  Affine2 res17(0.0);
  Affine2 res18(0.0);
  Affine2 res19(0.0);
  Affine2 res20(0.0);
  Affine2 res21(0.0);
  Affine2 res22(0.0);
  Affine2 res23(0.0);
  Affine2 res24(0.0);
  Affine2 res25(0.0);
  Affine2 res26(0.0);
  Affine2 res27(0.0);
  Affine2 res28(0.0);
  Affine2 res29(0.0);
  Affine2 res30(0.0);
  Affine2 res31(0.0);
  Affine2 res32(0.0);
  Affine2 res33(0.0);
  Affine2 res34(0.0);
  Affine2 res35(0.0);
  Affine2 res36(0.0);
  Affine2 res37(0.0);
  Affine2 res38(0.0);
  Affine2 res39(0.0);
  Affine2 res40(0.0);
  Affine2 res41(0.0);
  Affine2 res42(0.0);
  Affine2 res43(0.0);
  Affine2 res44(0.0);
  Affine2 res45(0.0);
  Affine2 res46(0.0);
  Affine2 res47(0.0);
  Affine2 res48(0.0);
  std::vector<int> key0;
  std::vector<int> key1;
  std::vector<int> key2;
  std::vector<int> key3;
  std::vector<int> key4;
  std::vector<int> key5;
  std::vector<int> key6;
  for (int j1 = 0; j1 < nbvar; ++j1) {
    // Common Factor

    key0.clear();
    key0.push_back (j1);
    Affine2 temp6 = edfr->eval_frechet(0, key0, y);

    key0.clear();
    key0.push_back (j); key0.push_back (j1);
    Affine2 temp35 = edfr->eval_frechet(1, key0, y);

    for (int j2 = 0; j2 < nbvar; ++j2) {
      // Common Factor

      key0.clear();
      key0.push_back (j2);
      Affine2 temp5 = edfr->eval_frechet(0, key0, y);

      key0.clear();
      key0.push_back (j); key0.push_back (j1); key0.push_back (j2);
      Affine2 temp22 = edfr->eval_frechet(2, key0, y);

      key0.clear();
      key0.push_back (j1); key0.push_back (j2);
      Affine2 temp39 = edfr->eval_frechet(1, key0, y);

      for (int j3 = 0; j3 < nbvar; ++j3) {
        // Common Factor

        key0.clear();
        key0.push_back (j3);
        Affine2 temp4 = edfr->eval_frechet(0, key0, y);

        key0.clear();
        key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp16 = edfr->eval_frechet(3, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j2); key0.push_back (j3);
        Affine2 temp38 = edfr->eval_frechet(2, key0, y);

        key0.clear();
        key0.push_back (j1); key0.push_back (j3);
        Affine2 temp30 = edfr->eval_frechet(1, key0, y);

        key0.clear();
        key0.push_back (j2); key0.push_back (j3);
        Affine2 temp25 = edfr->eval_frechet(1, key0, y);

        for (int j4 = 0; j4 < nbvar; ++j4) {
          // Common Factor

          key0.clear();
          key0.push_back (j4);
          Affine2 temp3 = edfr->eval_frechet(0, key0, y);

          key0.clear();
          key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp11 = edfr->eval_frechet(4, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp37 = edfr->eval_frechet(3, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp33 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j1); key0.push_back (j4);
          Affine2 temp28 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j3); key0.push_back (j4);
          Affine2 temp24 = edfr->eval_frechet(2, key0, y);

          key0.clear();
          key0.push_back (j2); key0.push_back (j4);
          Affine2 temp20 = edfr->eval_frechet(1, key0, y);

          key0.clear();
          key0.push_back (j3); key0.push_back (j4);
          Affine2 temp18 = edfr->eval_frechet(1, key0, y);

          for (int j5 = 0; j5 < nbvar; ++j5) {
            // Common Factor

            key0.clear();
            key0.push_back (j5);
            Affine2 temp2 = edfr->eval_frechet(0, key0, y);

            key0.clear();
            key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp9 = edfr->eval_frechet(5, key0, y);

            key0.clear();
            key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp36 = edfr->eval_frechet(4, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp23 = edfr->eval_frechet(3, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp31 = edfr->eval_frechet(2, key0, y);

            key0.clear();
            key0.push_back (j2); key0.push_back (j5);
            Affine2 temp27 = edfr->eval_frechet(1, key0, y);

            key0.clear();
            key0.push_back (j3); key0.push_back (j4); key0.push_back (j5);
            Affine2 temp17 = edfr->eval_frechet(2, key0, y);

            key0.clear();
            key0.push_back (j3); key0.push_back (j5);
            Affine2 temp14 = edfr->eval_frechet(1, key0, y);

            key0.clear();
            key0.push_back (j4); key0.push_back (j5);
            Affine2 temp12 = edfr->eval_frechet(1, key0, y);

            for (int j6 = 0; j6 < nbvar; ++j6) {
              // Common Factor

              key0.clear();
              key0.push_back (j6);
              Affine2 temp1 = edfr->eval_frechet(0, key0, y);

              key0.clear();
              key0.push_back (j); key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp7 = edfr->eval_frechet(6, key0, y);

              key0.clear();
              key0.push_back (j1); key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp34 = edfr->eval_frechet(5, key0, y);

              key0.clear();
              key0.push_back (j2); key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp21 = edfr->eval_frechet(4, key0, y);

              key0.clear();
              key0.push_back (j2); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp29 = edfr->eval_frechet(3, key0, y);

              key0.clear();
              key0.push_back (j2); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp32 = edfr->eval_frechet(2, key0, y);

              key0.clear();
              key0.push_back (j3); key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp15 = edfr->eval_frechet(3, key0, y);

              key0.clear();
              key0.push_back (j3); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp19 = edfr->eval_frechet(2, key0, y);

              key0.clear();
              key0.push_back (j3); key0.push_back (j6);
              Affine2 temp26 = edfr->eval_frechet(1, key0, y);

              key0.clear();
              key0.push_back (j4); key0.push_back (j5); key0.push_back (j6);
              Affine2 temp10 = edfr->eval_frechet(2, key0, y);

              key0.clear();
              key0.push_back (j4); key0.push_back (j6);
              Affine2 temp13 = edfr->eval_frechet(1, key0, y);

              key0.clear();
              key0.push_back (j5); key0.push_back (j6);
              Affine2 temp8 = edfr->eval_frechet(1, key0, y);

              // Tree computation
              res1 += temp1 * temp2 * temp3 * temp4 * temp5 * temp6 * temp7;
              res2 += temp1 * temp8 * temp3 * temp4 * temp5 * temp6 * temp9;
              res3 += temp1 * temp2 * temp10 * temp4 * temp5 * temp6 * temp11;
              res4 += temp1 * temp8 * temp12 * temp4 * temp5 * temp6 * temp11;
              res5 += temp1 * temp13 * temp2 * temp14 * temp5 * temp6 * temp11;
              res6 += temp1 * temp2 * temp3 * temp15 * temp5 * temp6 * temp16;
              res7 += temp1 * temp8 * temp3 * temp17 * temp5 * temp6 * temp16;
              res8 += temp1 * temp2 * temp10 * temp18 * temp5 * temp6 * temp16;
              res9 += temp1 * temp8 * temp12 * temp18 * temp5 * temp6 * temp16;
              res10 += temp1 * temp2 * temp19 * temp3 * temp20 * temp6 * temp16;
              res11 += temp1 * temp8 * temp14 * temp3 * temp20 * temp6 * temp16;
              res12 += temp1 * temp2 * temp3 * temp4 * temp21 * temp6 * temp22;
              res13 += temp1 * temp8 * temp3 * temp4 * temp23 * temp6 * temp22;
              res14 += temp1 * temp2 * temp10 * temp4 * temp24 * temp6 * temp22;
              res15 += temp1 * temp8 * temp12 * temp4 * temp24 * temp6 * temp22;
              res16 += temp1 * temp13 * temp2 * temp14 * temp24 * temp6 * temp22;
              res17 += temp1 * temp2 * temp3 * temp15 * temp25 * temp6 * temp22;
              res18 += temp1 * temp8 * temp3 * temp17 * temp25 * temp6 * temp22;
              res19 += temp1 * temp2 * temp10 * temp18 * temp25 * temp6 * temp22;
              res20 += temp1 * temp8 * temp12 * temp18 * temp25 * temp6 * temp22;
              res21 += temp1 * temp26 * temp2 * temp27 * temp3 * temp28 * temp16;
              res22 += temp1 * temp2 * temp3 * temp29 * temp4 * temp30 * temp22;
              res23 += temp1 * temp8 * temp3 * temp31 * temp4 * temp30 * temp22;
              res24 += temp1 * temp2 * temp10 * temp20 * temp4 * temp30 * temp22;
              res25 += temp1 * temp8 * temp12 * temp20 * temp4 * temp30 * temp22;
              res26 += temp1 * temp2 * temp32 * temp3 * temp4 * temp33 * temp22;
              res27 += temp1 * temp8 * temp27 * temp3 * temp4 * temp33 * temp22;
              res28 += temp1 * temp2 * temp3 * temp4 * temp5 * temp34 * temp35;
              res29 += temp1 * temp8 * temp3 * temp4 * temp5 * temp36 * temp35;
              res30 += temp1 * temp2 * temp10 * temp4 * temp5 * temp37 * temp35;
              res31 += temp1 * temp8 * temp12 * temp4 * temp5 * temp37 * temp35;
              res32 += temp1 * temp13 * temp2 * temp14 * temp5 * temp37 * temp35;
              res33 += temp1 * temp2 * temp3 * temp15 * temp5 * temp38 * temp35;
              res34 += temp1 * temp8 * temp3 * temp17 * temp5 * temp38 * temp35;
              res35 += temp1 * temp2 * temp10 * temp18 * temp5 * temp38 * temp35;
              res36 += temp1 * temp8 * temp12 * temp18 * temp5 * temp38 * temp35;
              res37 += temp1 * temp8 * temp27 * temp3 * temp18 * temp30 * temp22;
              res38 += temp1 * temp2 * temp19 * temp3 * temp20 * temp38 * temp35;
              res39 += temp1 * temp8 * temp14 * temp3 * temp20 * temp38 * temp35;
              res40 += temp1 * temp2 * temp3 * temp4 * temp21 * temp39 * temp35;
              res41 += temp1 * temp8 * temp3 * temp4 * temp23 * temp39 * temp35;
              res42 += temp1 * temp2 * temp10 * temp4 * temp24 * temp39 * temp35;
              res43 += temp1 * temp8 * temp12 * temp4 * temp24 * temp39 * temp35;
              res44 += temp1 * temp13 * temp2 * temp14 * temp24 * temp39 * temp35;
              res45 += temp1 * temp2 * temp3 * temp15 * temp25 * temp39 * temp35;
              res46 += temp1 * temp8 * temp3 * temp17 * temp25 * temp39 * temp35;
              res47 += temp1 * temp2 * temp10 * temp18 * temp25 * temp39 * temp35;
              res48 += temp1 * temp8 * temp12 * temp18 * temp25 * temp39 * temp35;


	      res1.compact();
	      res2.compact();
	      res3.compact();
	      res4.compact();
	      res5.compact();
	      res6.compact();
	      res7.compact();
	      res8.compact();
	      res9.compact();
	      res10.compact();
	      res11.compact();
	      res12.compact();
	      res13.compact();
	      res14.compact();
	      res15.compact();
	      res16.compact();
	      res17.compact();
	      res18.compact();
	      res19.compact();
	      res20.compact();
	      res21.compact();
	      res22.compact();
	      res23.compact();
	      res24.compact();
	      res25.compact();
	      res26.compact();
	      res27.compact();
	      res28.compact();
	      res29.compact();
	      res30.compact();
	      res31.compact();
	      res32.compact();
	      res33.compact();
	      res34.compact();
	      res35.compact();
	      res36.compact();
	      res37.compact();
	      res38.compact();
	      res39.compact();
	      res40.compact();
	      res41.compact();
	      res42.compact();
	      res43.compact();
	      res44.compact();
	      res45.compact();
	      res46.compact();
	      res47.compact();
	      res48.compact();

            }
          }
        }
      }
    }
  }

return Interval(0.24999999999977263E-2, 0.25000000000020561E-2) * res1 + Interval(0.37499999999965894E-1, 0.37500000000030842E-1) * res2 + Interval(0.49999999999945643E-1, 0.50000000000041123E-1) * res3 + Interval(0.49999999999941203E-1, 0.50000000000056666E-1) * res4 + Interval(0.11249999999988769E0, 0.11250000000010252E0) * res5 + Interval(-0.50000000000044453E-1, -0.49999999999964517E-1) * res6 + Interval(-0.15000000000013336E0, -0.14999999999987357E0) * res7 + Interval(-0.50000000000054445E-1, -0.49999999999947864E-1) * res8 + Interval(-0.50000000000064437E-1, -0.49999999999924549E-1) * res9 + Interval(0.14999999999983027E0, 0.15000000000013669E0) * res10 + Interval(0.14999999999982361E0, 0.15000000000017666E0) * res11 + Interval(-0.37500000000017852E-1, -0.37499999999983213E-1) * res12 + Interval(-0.22500000000011511E0, -0.22499999999989928E0) * res13 + Interval(-0.15000000000008740E0, -0.14999999999991687E0) * res14 + Interval(-0.15000000000010338E0, -0.14999999999989022E0) * res15 + Interval(-0.11250000000006555E0, -0.11249999999994165E0) * res16 + Interval(0.49999999999968514E-1, 0.50000000000031131E-1) * res17 + Interval(0.14999999999989755E0, 0.15000000000010538E0) * res18 + Interval(0.49999999999952527E-1, 0.50000000000051781E-1) * res19 + Interval(0.49999999999938538E-1, 0.50000000000071099E-1) * res20 + Interval(0.37499999999959233E-1, 0.37500000000039169E-1) * res21 + Interval(-0.50000000000044453E-1, -0.49999999999961187E-1) * res22 + Interval(-0.15000000000013336E0, -0.14999999999986358E0) * res23 + Interval(-0.50000000000057776E-1, -0.49999999999944533E-1) * res24 + Interval(-0.50000000000071099E-1, -0.49999999999921219E-1) * res25 + Interval(0.24999999999972822E-1, 0.25000000000023892E-1) * res26 + Interval(0.49999999999936762E-1, 0.50000000000056666E-1) * res27 + Interval(-0.15000000000004121E-1, -0.14999999999996128E-1) * res28 + Interval(-0.15000000000004565E0, -0.14999999999995905E0) * res29 + Interval(-0.15000000000004565E0, -0.14999999999995461E0) * res30 + Interval(-0.15000000000005009E0, -0.14999999999994351E0) * res31 + Interval(-0.22500000000006848E0, -0.22499999999993192E0) * res32 + Interval(0.99999999999973443E-1, 0.10000000000002618E0) * res33 + Interval(0.29999999999990867E0, 0.30000000000009019E0) * res34 + Interval(0.99999999999961786E-1, 0.10000000000004228E0) * res35 + Interval(0.99999999999950129E-1, 0.10000000000005671E0) * res36 + Interval(0.24999999999966160E-1, 0.25000000000034994E-1) * res37 + Interval(-0.15000000000004787E0, -0.14999999999995239E0) * res38 + Interval(-0.15000000000005453E0, -0.14999999999993907E0) * res39 + Interval(0.37499999999991207E-1, 0.37500000000008638E-1) * res40 + Interval(0.22499999999994391E0, 0.22500000000005849E0) * res41 + Interval(0.14999999999995861E0, 0.15000000000004521E0) * res42 + Interval(0.14999999999994929E0, 0.15000000000005542E0) * res43 + Interval(0.11249999999997029E0, 0.11250000000003124E0) * res44 + Interval(-0.50000000000016698E-1, -0.49999999999981615E-1) * res45 + Interval(-0.15000000000005542E0, -0.14999999999993752E0) * res46 + Interval(-0.50000000000026024E-1, -0.49999999999970512E-1) * res47 + Interval(-0.50000000000034239E-1, -0.49999999999961631E-1) * res48;

}

}
