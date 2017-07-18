program dbg
#ifdef __LCTERM__
write (*,*) "LC is defined"
#else
write (*,*) "LC is NOT defined"
#endif
stop
end program dbg
