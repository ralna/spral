/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * Licence: BSD licence, see LICENCE file for details
 *
 */
#pragma once

#include <cstdio>
#include <iostream>
#include <ios>
#include <stdexcept>

#define ANSI_COLOR_RED     "\x1b[31;1m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_BLUE    "\x1b[34;1m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define TEST(func) \
if(!func) { \
   printf(ANSI_COLOR_BLUE "%-20s" ANSI_COLOR_RESET " [ " ANSI_COLOR_GREEN "ok" ANSI_COLOR_RESET " ]\n", #func); \
} else { \
   printf(ANSI_COLOR_BLUE "%-20s" ANSI_COLOR_RESET " [" ANSI_COLOR_RED "fail" ANSI_COLOR_RESET "]\n", #func); \
   nerr++; \
}

namespace spral { namespace test {

class AssertError : public std::runtime_error {
public:
   AssertError()
   : runtime_error("Assertion Failed")
   {}
};

class Expectation {
public:
   Expectation(bool test)
   : passed_(test)
   {}
   template<typename T>
   Expectation const& operator<< (T in) const {
      if(!passed_) std::cout << in; // Only print stuff if we've failed the test
      return *this;
   }
   Expectation const& operator<< (std::ostream& (*in)(std::ostream&)) const {
      if(!passed_) std::cout << in; // Only print stuff if we've failed the test
      return *this;
   }
   /*Expectation const& operator<< (std::basic_ios& (*in)(std::basic_ios&)) const {
      if(!passed_) std::cout << in; // Only print stuff if we've failed the test
      return *this;
   }*/
   Expectation const& operator<< (std::ios_base& (*in)(std::ios_base&)) const {
      if(!passed_) std::cout << in; // Only print stuff if we've failed the test
      return *this;
   }
private:
   bool passed_;
};

class Assertion : public Expectation {
public:
   Assertion(bool test)
   : Expectation(test) {
      throw AssertError();
   }
};

}} /* namespace spral::test */

#define ASSERT_TRUE(a) \
if(!(a)) { \
   printf("Test failure at %s:%d\n", __FILE__, __LINE__); \
   printf("ASSERT_TRUE(%s) failed\n", #a); \
   return 1; \
}

#define ASSERT_EQ(a, b) \
if((a) != (b)) { \
   printf("Test failure at %s:%d\n", __FILE__, __LINE__); \
   printf("ASSERT_EQ(%s, %s) failed\n", #a, #b); \
   return 1; \
}

#define ASSERT_LE(a, b) \
if((a) > (b)) { \
   printf("Test failure at %s:%d\n", __FILE__, __LINE__); \
   printf("ASSERT_LE(%s, %s)\n", #a, #b); \
   return 1; \
}

#define EXPECT_LE(a, b) \
   if(!((a) <= (b))) failed=true; \
   spral::test::Expectation((a) <= (b)) \
   << "Test failure at " __FILE__ ":" << __LINE__ << std::endl \
   << "EXPECT_LE(" #a ", " #b ") : " << (a) << " > " << (b) << std::endl

#define EXPECT_EQ(a, b) \
   if(!((a) == (b))) failed=true; \
   spral::test::Expectation((a) == (b)) \
   << "Test failure at " __FILE__ ":" << __LINE__ << std::endl \
   << "EXPECT_EQ(" #a ", " #b ") : " << (a) << " != " << (b)  << std::endl

void gen_posdef(int n, double* a, int lda);
void gen_sym_indef(int n, double* a, int lda);
void gen_rhs(int n, double const* a, int lda, double* rhs);
void print_vec(char const* format, int n, double const* vec);
void print_mat(char const* format, int n, double const* a, int lda, int *perm=nullptr);
double backward_error(int n, double const* a, int lda, double const* rhs, int nrhs, double const* soln, int ldsoln);
double forward_error(int n, int nrhs, double const* soln, int ldx);
