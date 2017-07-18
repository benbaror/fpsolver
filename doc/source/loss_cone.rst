**************
Loss Cone term
**************



Bahcall and Wolf Loss Cone term
===============================
Bahcall and Wolf use a loss cone term based on Lighman and Shapiro
(1977) analysis.

.. math::

   R_{NR} =
   \frac{0.5g(x,\tau)x^{5/2}}{\alpha\ln(Q))[\beta+q(x,\tau)^{-1}\ln(x_D/4x)]}

where  :math:`\alpha = Q^{-2}n_hr_h^3x_D`,

       :math:`Q = M_\bullet/M_\star`,

       :math:`q(x,\tau) = 1.6\ln(6Q/\pi)\alpha x^{-5/2}g(x,\tau)`,

       :math:`\beta = \frac{1}{2}\int_0^\infty dy 
       \frac{e^{-y}}{e^{-y} + \sqrt{\pi y} [1+erf(\sqrt{y})]} \approx
       1.8`

The same term has been used by Hopman and Alexander 2006A.


In the empty lose cone regime we have

.. math::

   R_{NR} =
   \frac{0.8g(x,\tau)^2\ln(6Q/\pi)}{\ln(Q)ln(x_D/4x)}
   \approx
   \frac{0.8g(x,\tau)^2}{ln(x_D/4x)}

Hopman and Alexander Loss Cone term
===================================

Hopman and Alexander 2006B and Alexander  and Hopman 2009 used a
simpler loss cone term for the NR.

.. math::

	R_{NR} = \frac{g_M(x)}{\tau_r(x)\ln[J_c(x)/J_{lc}]}

where :math:`\tau_r(x) = \frac{M^2_\star}{\sum_M g_M(x)M^2}`,
	      
      :math:`Jlso = 4GM_\bullet/c`

Resonant Relaxation Loss Cone term
==================================

Hopman and Alexander 2006B add a RR loss cone term

.. math::

	R_{RR} = g(x,\tau)/\tau_{RR}(x)

where :math:`\tau_{RR} \equiv \bar{T}_{RR}^S(x)/T_h`

Compering the different treatments
==================================

.. plot:: source/pyplots/loss_cone_comp.py


Compering the Fokker Plank results of the different treatments
==============================================================

.. plot:: source/pyplots/loss_cone_results.py
