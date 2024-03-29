function dTdLx = OligopolyODELx(Lx,T,G,K,Ly,M,px,py)
%OLIGOPOLYODELX
%    [DT(1)LX,DT(2)LX] = OLIGOPOLYODELX(G,K,LX,LY,M,T(1),T(2),PX,PY)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    30-Jun-2020 18:08:49

dTdLx = zeros(2,1);

t2 = Lx.*T(2);
t3 = Ly.*T(1);
t4 = Lx.*Ly;
t5 = t2+t3+t4;
t6 = 1.0./t5;
t7 = Ly.*T(1).*py.*t6;
t8 = py-1.0;
t9 = Lx.*T(2).*t6.*t8;
t10 = -py+t7+t9;
t11 = px-1.0;
t12 = M.^2;
t13 = t12-4.0;
t14 = 1.0./t13;
t15 = 1.0./t5.^2;
t16 = 1.0./t5.^3;
t17 = Lx.*Ly.*py.*t15;
t18 = Lx.*Ly.*t8.*t15;
t19 = Lx.^2;
t20 = Ly.^2;
t22 = Ly.*T(2).*t8.*t16.*t19.*2.0;
t23 = Lx.*T(1).*py.*t16.*t20.*2.0;
t21 = t17+t18-t22-t23;
t24 = Lx.*Ly.*px.*t15;
t25 = Lx.*Ly.*t11.*t15;
t27 = Lx.*T(1).*t11.*t16.*t20.*2.0;
t28 = Ly.*T(2).*px.*t16.*t19.*2.0;
t26 = t24+t25-t27-t28;
t29 = T(1).*py.*t15.*t20;
t30 = Lx.*Ly.*T(2).*t8.*t15;
t44 = Ly.*py.*t6;
t31 = t29+t30-t44;
t32 = T(2).*t8.*t15.*t19;
t33 = Lx.*Ly.*T(1).*py.*t15;
t35 = Lx.*t6.*t8;
t34 = t32+t33-t35;
t36 = T(2).*px.*t15.*t19;
t37 = Lx.*Ly.*T(1).*t11.*t15;
t39 = Lx.*px.*t6;
t38 = t36+t37-t39;
t40 = K.*t34.*2.0;
t41 = K.*M.*t38;
t42 = Lx.*px.*t6.*2.0;
t85 = M.*t34;
t86 = T(2).*px.*t15.*t19.*2.0;
t87 = Lx.*Ly.*T(1).*t11.*t15.*2.0;
t43 = t40+t41+t42-t85-t86-t87;
t45 = K.*t31.*2.0;
t46 = Ly.*t6.*t11.*2.0;
t47 = T(1).*t11.*t15.*t20;
t48 = Lx.*Ly.*T(2).*px.*t15;
t52 = Ly.*t6.*t11;
t49 = t47+t48-t52;
t50 = K.*M.*t49;
t94 = M.*t31;
t95 = T(1).*t11.*t15.*t20.*2.0;
t96 = Lx.*Ly.*T(2).*px.*t15.*2.0;
t51 = t45+t46+t50-t94-t95-t96;
t53 = px.*2.0;
t54 = K.*t10.*2.0;
t55 = Lx.*T(2).*px.*t6;
t56 = Ly.*T(1).*t6.*t11;
t57 = -px+t55+t56;
t58 = K.*M.*t57;
t100 = M.*t10;
t101 = Lx.*T(2).*px.*t6.*2.0;
t102 = Ly.*T(1).*t6.*t11.*2.0;
t59 = t53+t54+t58-t100-t101-t102;
t60 = K.*t21.*2.0;
t61 = K.*M.*t26;
t62 = Lx.*T(1).*t11.*t16.*t20.*4.0;
t63 = Ly.*T(2).*px.*t16.*t19.*4.0;
t74 = M.*t21;
t75 = Lx.*Ly.*px.*t15.*2.0;
t76 = Lx.*Ly.*t11.*t15.*2.0;
t64 = t60+t61+t62+t63-t74-t75-t76;
t65 = py.*2.0;
t66 = K.*t57.*2.0;
t67 = K.*M.*t10;
t97 = M.*t57;
t98 = Ly.*T(1).*py.*t6.*2.0;
t99 = Lx.*T(2).*t6.*t8.*2.0;
t68 = t65+t66+t67-t97-t98-t99;
t69 = K.*t26.*2.0;
t70 = K.*M.*t21;
t71 = Ly.*T(2).*t8.*t16.*t19.*4.0;
t72 = Lx.*T(1).*py.*t16.*t20.*4.0;
t103 = M.*t26;
t104 = Lx.*Ly.*py.*t15.*2.0;
t105 = Lx.*Ly.*t8.*t15.*2.0;
t73 = t69+t70+t71+t72-t103-t104-t105;
t77 = K.*t49.*2.0;
t78 = K.*M.*t31;
t79 = Ly.*py.*t6.*2.0;
t91 = M.*t49;
t92 = T(1).*py.*t15.*t20.*2.0;
t93 = Lx.*Ly.*T(2).*t8.*t15.*2.0;
t80 = t77+t78+t79-t91-t92-t93;
t81 = K.*t38.*2.0;
t82 = Lx.*t6.*t8.*2.0;
t83 = K.*M.*t34;
t88 = M.*t38;
t89 = T(2).*t8.*t15.*t19.*2.0;
t90 = Lx.*Ly.*T(1).*py.*t15.*2.0;
t84 = t81+t82+t83-t88-t89-t90;
t106 = Lx.*T(2).*px.*t16.*t19.*2.0;
t107 = Ly.*T(1).*t11.*t16.*t19.*2.0;
t112 = px.*t15.*t19.*2.0;
t108 = t106+t107-t112;
t109 = K.*t57;
t110 = t14.*t68;
t147 = M.*t14.*t59;
t111 = py-t7-t9+t109+t110-t147;
t113 = K.*t108.*2.0;
t114 = t8.*t15.*t19.*4.0;
t115 = Lx.*T(2).*t8.*t16.*t19.*2.0;
t116 = Ly.*T(1).*py.*t16.*t19.*2.0;
t120 = t8.*t15.*t19.*2.0;
t117 = t115+t116-t120;
t118 = K.*M.*t117;
t189 = M.*t108;
t190 = Lx.*T(2).*t8.*t16.*t19.*4.0;
t191 = Ly.*T(1).*py.*t16.*t19.*4.0;
t119 = t113+t114+t118-t189-t190-t191;
t121 = M.*t14.*t43;
t141 = t14.*t84;
t142 = K.*t38;
t122 = t32+t33-t35+t121-t141-t142;
t123 = G.*2.0;
t124 = Ly.*T(1).*py.*t16.*t20.*2.0;
t125 = Lx.*T(2).*t8.*t16.*t20.*2.0;
t130 = py.*t15.*t20.*2.0;
t126 = t124+t125-t130;
t127 = K.*t10;
t128 = t14.*t59;
t185 = M.*t14.*t68;
t129 = px-t55-t56+t127+t128-t185;
t131 = K.*t126.*2.0;
t132 = t11.*t15.*t20.*4.0;
t133 = Ly.*T(1).*t11.*t16.*t20.*2.0;
t134 = Lx.*T(2).*px.*t16.*t20.*2.0;
t138 = t11.*t15.*t20.*2.0;
t135 = t133+t134-t138;
t136 = K.*M.*t135;
t242 = M.*t126;
t243 = Ly.*T(1).*t11.*t16.*t20.*4.0;
t244 = Lx.*T(2).*px.*t16.*t20.*4.0;
t137 = t131+t132+t136-t242-t243-t244;
t139 = M.*t14.*t80;
t183 = t14.*t51;
t184 = K.*t31;
t140 = t47+t48-t52+t139-t183-t184;
t143 = Ly+T(2);
t144 = Lx.*T(2).*px.*t15.*t143;
t145 = Ly.*T(1).*t11.*t15.*t143;
t152 = T(2).*px.*t6;
t146 = t144+t145-t152;
t148 = px.*t6;
t149 = T(2).*px.*t16.*t19.*t143.*2.0;
t150 = Lx.*Ly.*T(1).*t11.*t16.*t143.*2.0;
t161 = Lx.*T(2).*px.*t15.*2.0;
t162 = Ly.*T(1).*t11.*t15;
t163 = Lx.*px.*t15.*t143;
t151 = t148+t149+t150-t161-t162-t163;
t153 = K.*t146.*2.0;
t154 = T(2).*t6.*t8.*2.0;
t155 = Ly.*T(1).*py.*t15.*t143;
t156 = Lx.*T(2).*t8.*t15.*t143;
t160 = T(2).*t6.*t8;
t157 = t155+t156-t160;
t158 = K.*M.*t157;
t210 = M.*t146;
t211 = Ly.*T(1).*py.*t15.*t143.*2.0;
t212 = Lx.*T(2).*t8.*t15.*t143.*2.0;
t159 = t153+t154+t158-t210-t211-t212;
t164 = K.*t151.*2.0;
t165 = t6.*t8;
t166 = T(2).*t8.*t16.*t19.*t143.*2.0;
t167 = Lx.*Ly.*T(1).*py.*t16.*t143.*2.0;
t174 = Ly.*T(1).*py.*t15;
t175 = Lx.*t8.*t15.*t143;
t176 = Lx.*T(2).*t8.*t15.*2.0;
t168 = t165+t166+t167-t174-t175-t176;
t169 = K.*M.*t168;
t170 = Ly.*T(1).*py.*t15.*2.0;
t171 = Lx.*t8.*t15.*t143.*2.0;
t172 = Lx.*T(2).*t8.*t15.*4.0;
t283 = M.*t151;
t284 = t6.*t8.*2.0;
t285 = T(2).*t8.*t16.*t19.*t143.*4.0;
t286 = Lx.*Ly.*T(1).*py.*t16.*t143.*4.0;
t173 = t164+t169+t170+t171+t172-t283-t284-t285-t286;
t177 = t14.*t64;
t178 = K.*t21;
t225 = M.*t14.*t73;
t179 = -t24-t25+t27+t28+t177+t178-t225;
t180 = t14.*t59.*t179;
t181 = M.*t14.*t84;
t226 = t14.*t43;
t227 = K.*t34;
t182 = t14.*t51.*(t36+t37-t39+t181-t226-t227);
t186 = t14.*t64.*t129;
t187 = t14.*t43.*t140;
t188 = t180+t182+t186+t187;
t192 = t14.*t111.*t119;
t193 = K.*t117.*2.0;
t194 = px.*t15.*t19.*4.0;
t195 = K.*M.*t108;
t237 = M.*t117;
t238 = Lx.*T(2).*px.*t16.*t19.*4.0;
t239 = Ly.*T(1).*t11.*t16.*t19.*4.0;
t196 = t193+t194+t195-t237-t238-t239;
t197 = M.*t14.*t196;
t235 = t14.*t119;
t236 = K.*t108;
t198 = t115+t116-t120+t197-t235-t236;
t240 = t14.*t68.*t198;
t241 = t14.*t84.*t122.*2.0;
t199 = t123+t192-t240-t241;
t200 = Ly.*T(2).*t8.*t15;
t201 = Ly.*py.*t15.*t143;
t213 = T(1).*py.*t16.*t20.*t143.*2.0;
t214 = Lx.*Ly.*T(2).*t8.*t16.*t143.*2.0;
t202 = t200+t201-t213-t214;
t203 = K.*t157.*2.0;
t204 = T(2).*px.*t6.*2.0;
t205 = K.*M.*t146;
t207 = M.*t157;
t208 = Lx.*T(2).*px.*t15.*t143.*2.0;
t209 = Ly.*T(1).*t11.*t15.*t143.*2.0;
t206 = t203+t204+t205-t207-t208-t209;
t215 = K.*t202.*2.0;
t216 = Ly.*T(2).*px.*t15;
t217 = Ly.*t11.*t15.*t143;
t223 = T(1).*t11.*t16.*t20.*t143.*2.0;
t224 = Lx.*Ly.*T(2).*px.*t16.*t143.*2.0;
t218 = t216+t217-t223-t224;
t219 = K.*M.*t218;
t220 = T(1).*t11.*t16.*t20.*t143.*4.0;
t221 = Lx.*Ly.*T(2).*px.*t16.*t143.*4.0;
t258 = M.*t202;
t259 = Ly.*T(2).*px.*t15.*2.0;
t260 = Ly.*t11.*t15.*t143.*2.0;
t222 = t215+t219+t220+t221-t258-t259-t260;
t228 = t14.*t73;
t229 = K.*t26;
t273 = M.*t14.*t64;
t230 = -t17-t18+t22+t23+t228+t229-t273;
t231 = t14.*t68.*t230;
t232 = M.*t14.*t51;
t274 = t14.*t80;
t275 = K.*t49;
t233 = t14.*t84.*(t29+t30-t44+t232-t274-t275);
t234 = t14.*t73.*t111;
t245 = t14.*t129.*t137;
t246 = K.*t135.*2.0;
t247 = py.*t15.*t20.*4.0;
t248 = K.*M.*t126;
t278 = M.*t135;
t279 = Ly.*T(1).*py.*t16.*t20.*4.0;
t280 = Lx.*T(2).*t8.*t16.*t20.*4.0;
t249 = t246+t247+t248-t278-t279-t280;
t250 = M.*t14.*t249;
t276 = t14.*t137;
t277 = K.*t126;
t251 = t133+t134-t138+t250-t276-t277;
t281 = t14.*t59.*t251;
t282 = t14.*t51.*t140.*2.0;
t252 = t123+t245-t281-t282;
t253 = t199.*t252;
t254 = t14.*t80.*t122;
t255 = t231+t233+t234+t254;
t256 = t253-t188.*t255;
t257 = 1.0./t256;
t261 = t14.*t129.*t222;
t262 = M.*t14.*t159;
t263 = t14.*t51.*(t144+t145-t152+t262-K.*t157-t14.*t206);
t264 = t14.*t222;
t265 = K.*t202;
t266 = K.*t218.*2.0;
t267 = K.*M.*t202;
t268 = T(1).*py.*t16.*t20.*t143.*4.0;
t269 = Lx.*Ly.*T(2).*t8.*t16.*t143.*4.0;
t270 = t266+t267+t268+t269-M.*t218-Ly.*py.*t15.*t143.*2.0-Ly.*T(2).*t8.*t15.*2.0;
t271 = -t216-t217+t223+t224+t264+t265-M.*t14.*t270;
t272 = t14.*t59.*t271;
t287 = t14.*t111.*t173;
t288 = M.*t14.*t206;
t289 = t155+t156-t160+t288-K.*t146-t14.*t159;
t290 = K.*t151;
t291 = t14.*t173;
t292 = K.*t168.*2.0;
t293 = K.*M.*t151;
t294 = Lx.*T(2).*px.*t15.*4.0;
t295 = Ly.*T(1).*t11.*t15.*2.0;
t296 = Lx.*px.*t15.*t143.*2.0;
t297 = t292+t293+t294+t295+t296-M.*t168-px.*t6.*2.0-T(2).*px.*t16.*t19.*t143.*4.0-Lx.*Ly.*T(1).*t11.*t16.*t143.*4.0;
t298 = -t165-t166-t167+t174+t175+t176+t290+t291-M.*t14.*t297;
t299 = t14.*t68.*t298;
t300 = t287+t299-t14.*t122.*t159-t14.*t84.*t289;
t301 = t14.*t140.*t206;
t302 = t261+t263+t272+t301;
DTxLx = -t257.*(t188.*t300-t199.*t302);
DTyLx = -t257.*(t252.*t300-t255.*t302);

dTdLx(1) = DTxLx;
dTdLx(2) = DTyLx;
