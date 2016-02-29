/**
 * ThetaTaylor.java
 *
 * AMH11: Java implementation of the matrix exponential method
 *     described by Al-Mohy and Higham (2011)
 *
 * Copyright (c) 2016 Arman Bilge
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.compevol.amh11;

/**
 * @author Arman Bilge
 */
public final class ThetaTaylor {

    private ThetaTaylor() {}

    public static final double[] THETA = {
        0x1.0p-52,
        0x1.bb67ae3d95a85p-26,
        0x1.d12e5ffa932eap-17,
        0x1.64380e81b76e7p-12,
        0x1.3ab00ae0b2f86p-9,
        0x1.29103ce4c952p-7,
        0x1.86ab5053bf46p-6,
        0x1.98e1a7f5da6c9p-5,
        0x1.6ee8ec72ebb51p-4,
        0x1.2749677a06ce3p-3,
        0x1.b6c1434041e9p-3,
        0x1.32ce821b60379p-2,
        0x1.995f48227e91fp-2,
        0x1.071fd3a3ff9cfp-1,
        0x1.483c197d30f8ep-1,
        0x1.8f81d55d66fc9p-1,
        0x1.dc6ecd1b85f4cp-1,
        0x1.1742d83ebe5e8p0,
        0x1.42a8554bf53aap0,
        0x1.7031527aa9654p0,
        0x1.9fabd93841a48p0,
        0x1.d0ea79cba2882p0,
        0x1.01e20af98722bp1,
        0x1.1c09cb19130d1p1,
        0x1.36dbcacb0f3d6p1,
        0x1.524905a7f453ap1,
        0x1.6e43fb4a730d1p1,
        0x1.8ac087a4c98d6p1,
        0x1.a7b3be5db7d4p1,
        0x1.c513c9734d8d6p1,
        0x1.e2d7cb244a608p1,
        0x1.007be17e6659ap2,
        0x1.0fb63ae7c8696p2,
        0x1.1f17ac2e04b89p2,
        0x1.2e9d3e2f93f22p2,
        0x1.3e443fe710c35p2,
        0x1.4e0a3eff16a63p2,
        0x1.5ded014081364p2,
        0x1.6dea7ec18fd73p2,
        0x1.7e00dcbf13befp2,
        0x1.8e2e690b5465dp2,
        0x1.9e719600af0fbp2,
        0x1.aec8f6e819c38p2,
        0x1.bf333cc59e7dbp2,
        0x1.cfaf337e99f03p2,
        0x1.e03bbf4f0cecep2,
        0x1.f0d7da84a4c3fp2,
        0x1.00c149bba4f76p3,
        0x1.091d855c04f37p3,
        0x1.118038b80aa21p3,
        0x1.19e903f529223p3,
        0x1.22578e2158a21p3,
        0x1.2acb849b83e41p3,
        0x1.33449a8add296p3,
        0x1.3bc2886377fc8p3,
        0x1.44450b76b75bdp3,
        0x1.4ccbe58e4b696p3,
        0x1.5556dc909fe58p3,
        0x1.5de5ba2dbe0dap3,
        0x1.66784b93c1a27p3,
        0x1.6f0e612a1a66cp3,
        0x1.77a7ce52e9cf6p3,
        0x1.80446931e04afp3,
        0x1.88e40a780ed42p3,
        0x1.91868d3430bbdp3,
        0x1.9a2bcea6ff029p3,
        0x1.a2d3ae1b2a72ap3,
        0x1.ab7e0cc0a41efp3,
        0x1.b42acd8ae5303p3,
        0x1.bcd9d511ef1a1p3,
        0x1.c58b0975c4983p3,
        0x1.ce3e524422476p3,
        0x1.d6f39860436d8p3,
        0x1.dfaac5ec849dap3,
        0x1.e863c635ba755p3,
        0x1.f11e85a016bcdp3,
        0x1.f9daf19579db1p3,
        0x1.014c7c3a88e33p4,
        0x1.05ac44c21544bp4,
        0x1.0a0cca700b854p4,
        0x1.0e6e05b88d8efp4,
        0x1.12cfef6f2100ap4,
        0x1.173280c10724p4,
        0x1.1b95b32ff8f75p4,
        0x1.1ff9808d3f6b5p4,
        0x1.245de2f520a77p4,
        0x1.28c2d4ca9bc5dp4,
        0x1.2d2850b36d092p4,
        0x1.318e519455126p4,
        0x1.35f4d28d9e124p4,
        0x1.3a5bcef7da5e9p4,
        0x1.3ec34260d839dp4,
        0x1.432b2888c702p4,
        0x1.47937d5f8a8a1p4,
        0x1.4bfc3d023a8d6p4,
        0x1.506563b8cf8dfp4,
        0x1.54ceedf408e9p4,
        0x1.5938d84bb7a2ep4,
        0x1.5da31f7df3f48p4,
        0x1.620dc070b54ep4
    };

}
