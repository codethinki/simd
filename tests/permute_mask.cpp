/* SPDX-License-Identifier: BSD-3-Clause */
/* Copyright © 2026      GSI Helmholtzzentrum fuer Schwerionenforschung GmbH
 *                       Matthias Kretz <m.kretz@gsi.de>
 */

#include "unittest.h"
#include "complex_init.h"

template <typename V>
  struct Tests
  {
    using T = typename V::value_type;
    using M = typename V::mask_type;

    static constexpr V iota = test_iota<V, 1>;
    static constexpr T fill = T(-1);
    static constexpr M odd = M([](int i) { return 1 == (i & 1); });
    static constexpr M first = M([](int i) { return i == 0; });
    static constexpr M last = M([](int i) { return i == V::size() - 1; });
    static constexpr M lower_half = M([](int i) { return i < V::size() / 2; });

    ADD_TEST(compress_vec) {
      std::tuple {iota, M(true), M(false), odd, first, last},
      [](auto& t, V x, M all, M none, M k, M k0, M kn) {
	t.verify_equal(simd::compress(x, all), x);
	t.verify_equal(simd::compress(x, all, fill), x);
	t.verify_equal(simd::compress(x, none, fill), V(fill));
	t.verify_equal(simd::compress(x, k, fill),
		       V([](int i) { return 2 * i + 1 < V::size() ? iota[2 * i + 1] : fill; }));
	t.verify_equal(select(lower_half, simd::compress(x, k), V(fill)),
		       V([](int i) { return 2 * i + 1 < V::size() ? iota[2 * i + 1] : fill; }));
	t.verify_equal(simd::compress(x, k0, fill), V([](int i) { return i == 0 ? iota[0] : fill; }));
	t.verify_equal(simd::compress(x, kn, fill),
		       V([](int i) { return i == 0 ? iota[V::size() - 1] : fill; }));
      }
    };

    ADD_TEST(expand_vec) {
      std::tuple {iota, V(fill), M(true), M(false), odd},
      [](auto& t, V x, V orig, M all, M none, M k) {
	t.verify_equal(simd::expand(x, all), x);
	t.verify_equal(simd::expand(x, all, orig), x);
	t.verify_equal(simd::expand(x, none), V());
	t.verify_equal(simd::expand(x, none, orig), orig);
	t.verify_equal(simd::expand(x, k, orig),
		       V([](int i) { return 1 == (i & 1) ? iota[i / 2] : fill; }));
	t.verify_equal(simd::expand(simd::compress(x, k), k, orig), select(k, x, orig));
      }
    };

    ADD_TEST(compress_expand_mask) {
      std::tuple {odd, M(true), M(false), first},
      [](auto& t, M k, M all, M none, M k0) {
	t.verify_equal(simd::compress(k, all), k);
	t.verify_equal(simd::compress(k, none, true), M(true));
	t.verify_equal(simd::compress(k, k, false), lower_half);
	t.verify_equal(simd::compress(k0, k0, false), first);
	t.verify_equal(simd::compress(M(true), last, false), first);
	t.verify_equal(simd::expand(M(true), k), k);
	t.verify_equal(simd::expand(M(false), k, M(true)), !k);
	t.verify_equal(simd::expand(all, none, k), k);
      }
    };
  };
