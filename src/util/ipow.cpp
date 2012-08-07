int ipow(int base, int exp) {
	if (exp < 0) return 0;
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}
