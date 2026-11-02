/// Reports whether this C build contracts `a*b + c` into one fused
/// multiply-add. It is compiled with the same flags as the MAFFT sources in
/// this library, so FFI tests can tell whether the C side they compare
/// against follows the same contraction policy as the Rust build.
///
/// With a = 1 + 2^-30 and b = 1 - 2^-30, a*b = 1 - 2^-60 exactly. Rounded
/// to double it is 1.0, so the unfused sum with c = -1 is 0; the fused
/// result is -2^-60. The operands are read through `volatile` so the
/// expression is not folded at compile time.

int mafft_c_contracts_fma( void )
{
	volatile double va = 1.0 + 0x1p-30, vb = 1.0 - 0x1p-30, vc = -1.0;
	double a = va, b = vb, c = vc;
	double r = a * b + c;
	return r != 0.0;
}
