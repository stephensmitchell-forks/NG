//C++ TO C# CONVERTER TODO TASK: Check to ensure that the necessary preprocessor flags are defined:
internal static class DefineConstants
{
	public const double M_PI = 3.14159265358979323846;
	public const string PACKAGE_VERSION = "6.2-dev";
	public const int CLOCKS_PER_SEC = 1000000;
	public const int GZSTREAM_H = 1;
	public const double EPSGEOM = 1E-5;
	public const int ELEMENT_MAXPOINTS = 20;
	public const int ELEMENT2D_MAXPOINTS = 8;
	public const int MULTIPOINTGEOMINFO_MAX = 100;
	public const double DEFAULT_R = 0.0;
	public const double DEFAULT_G = 1.0;
	public const double DEFAULT_B = 0.0;
	public const int DEFAULT_BCNUM = 1;
	public const double DEFAULT_EPS = 2.5e-05;
#if OPENGLxx && WIN32
	public const int GL_CLAMP_TO_EDGE = 0x812F;
#endif
#if OPENGLxx && WIN32
	public const int GL_ARRAY_BUFFER = 0x8892;
#endif
#if OPENGLxx && WIN32
	public const int GL_ELEMENT_ARRAY_BUFFER = 0x8893;
#endif
#if OPENGLxx && WIN32
	public const int GL_STATIC_DRAW = 0x88E4;
#endif
#if OPENGLxx
	public const int STLBASE = 1;
#endif
#if WIN32
	public const char COMMASIGN = ':';
#else
	public const char COMMASIGN = ',';
#endif
}