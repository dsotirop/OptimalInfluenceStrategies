function dTdK = OligopolyODEK(K,T,G,Lx,Ly,M,px,py)
%OLIGOPOLYODEK
%    [DT(1)K,DT(2)K] = OLIGOPOLYODEK(G,K,LX,LY,M,T(1),T(2),PX,PY)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    30-Jun-2020 15:02:10

dTdK = zeros(2,1);

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
t152 = M.*t14.*t59;
t111 = py-t7-t9+t109+t110-t152;
t113 = K.*t108.*2.0;
t114 = t8.*t15.*t19.*4.0;
t115 = Lx.*T(2).*t8.*t16.*t19.*2.0;
t116 = Ly.*T(1).*py.*t16.*t19.*2.0;
t120 = t8.*t15.*t19.*2.0;
t117 = t115+t116-t120;
t118 = K.*M.*t117;
t156 = M.*t108;
t157 = Lx.*T(2).*t8.*t16.*t19.*4.0;
t158 = Ly.*T(1).*py.*t16.*t19.*4.0;
t119 = t113+t114+t118-t156-t157-t158;
t121 = M.*t14.*t43;
t154 = t14.*t84;
t155 = K.*t38;
t122 = t32+t33-t35+t121-t154-t155;
t123 = G.*2.0;
t124 = Ly.*T(1).*py.*t16.*t20.*2.0;
t125 = Lx.*T(2).*t8.*t16.*t20.*2.0;
t130 = py.*t15.*t20.*2.0;
t126 = t124+t125-t130;
t127 = K.*t10;
t128 = t14.*t59;
t149 = M.*t14.*t68;
t129 = px-t55-t56+t127+t128-t149;
t131 = K.*t126.*2.0;
t132 = t11.*t15.*t20.*4.0;
t133 = Ly.*T(1).*t11.*t16.*t20.*2.0;
t134 = Lx.*T(2).*px.*t16.*t20.*2.0;
t138 = t11.*t15.*t20.*2.0;
t135 = t133+t134-t138;
t136 = K.*M.*t135;
t197 = M.*t126;
t198 = Ly.*T(1).*t11.*t16.*t20.*4.0;
t199 = Lx.*T(2).*px.*t16.*t20.*4.0;
t137 = t131+t132+t136-t197-t198-t199;
t139 = M.*t14.*t80;
t147 = t14.*t51;
t148 = K.*t31;
t140 = t47+t48-t52+t139-t147-t148;
t141 = t14.*t64;
t142 = K.*t21;
t180 = M.*t14.*t73;
t143 = -t24-t25+t27+t28+t141+t142-t180;
t144 = t14.*t59.*t143;
t145 = M.*t14.*t84;
t181 = t14.*t43;
t182 = K.*t34;
t146 = t14.*t51.*(t36+t37-t39+t145-t181-t182);
t150 = t14.*t64.*t129;
t151 = -t42+t85+t86+t87;
t153 = -t53+t100+t101+t102;
t159 = t14.*t111.*t119;
t160 = K.*t117.*2.0;
t161 = px.*t15.*t19.*4.0;
t162 = K.*M.*t108;
t192 = M.*t117;
t193 = Lx.*T(2).*px.*t16.*t19.*4.0;
t194 = Ly.*T(1).*t11.*t16.*t19.*4.0;
t163 = t160+t161+t162-t192-t193-t194;
t164 = M.*t14.*t163;
t190 = t14.*t119;
t191 = K.*t108;
t165 = t115+t116-t120+t164-t190-t191;
t195 = t14.*t68.*t165;
t196 = t14.*t84.*t122.*2.0;
t166 = t123+t159-t195-t196;
t167 = T(1).*t15.*t20;
t168 = Lx.*Ly.*T(2).*t15;
t170 = Ly.*t6;
t169 = t167+t168-t170;
t171 = Lx.*T(2).*t6;
t172 = Ly.*T(1).*t6;
t173 = t171+t172-1.0;
t174 = K.*t169.*2.0;
t230 = M.*t169;
t175 = t174-t230;
t176 = K.*t173.*2.0;
t238 = M.*t173;
t177 = t176-t238;
t178 = t14.*t43.*t140;
t179 = t144+t146+t150+t178;
t183 = t14.*t73;
t184 = K.*t26;
t247 = M.*t14.*t64;
t185 = -t17-t18+t22+t23+t183+t184-t247;
t186 = t14.*t68.*t185;
t187 = M.*t14.*t51;
t248 = t14.*t80;
t249 = K.*t49;
t188 = t14.*t84.*(t29+t30-t44+t187-t248-t249);
t189 = t14.*t73.*t111;
t200 = t14.*t129.*t137;
t201 = K.*t135.*2.0;
t202 = py.*t15.*t20.*4.0;
t203 = K.*M.*t126;
t215 = M.*t135;
t216 = Ly.*T(1).*py.*t16.*t20.*4.0;
t217 = Lx.*T(2).*t8.*t16.*t20.*4.0;
t204 = t201+t202+t203-t215-t216-t217;
t205 = M.*t14.*t204;
t213 = t14.*t137;
t214 = K.*t126;
t206 = t133+t134-t138+t205-t213-t214;
t218 = t14.*t59.*t206;
t219 = t14.*t51.*t140.*2.0;
t207 = t123+t200-t218-t219;
t208 = t166.*t207;
t209 = t14.*t80.*t122;
t210 = t186+t188+t189+t209;
t211 = t208-t179.*t210;
t212 = 1.0./t211;
t220 = t14.*t151;
t221 = -t82+t88+t89+t90;
t222 = t36+t37-t39+t220-M.*t14.*t221;
t223 = t14.*t68.*t222;
t224 = t14.*t111.*t151;
t225 = t14.*t153;
t226 = -t65+t97+t98+t99;
t227 = -px+t55+t56+t225-M.*t14.*t226;
t228 = t14.*t84.*t227;
t229 = t223+t224+t228-t14.*t122.*t153;
t231 = t14.*t175;
t232 = K.*t169;
t233 = Ly.*t6.*2.0;
t234 = K.*M.*t169;
t235 = t233+t234-T(1).*t15.*t20.*2.0-Lx.*Ly.*T(2).*t15.*2.0;
t236 = t231+t232-M.*t14.*t235;
t237 = t14.*t59.*t236;
t239 = t14.*t177;
t240 = K.*t173;
t241 = K.*M.*t173;
t242 = t241-Lx.*T(2).*t6.*2.0-Ly.*T(1).*t6.*2.0+2.0;
t243 = t239+t240-M.*t14.*t242;
t244 = t14.*t51.*t243;
t245 = t14.*t129.*t175;
t246 = t237+t244+t245-t14.*t140.*t177;
DTxK = t212.*(t179.*t229+t166.*t246);
DTyK = t212.*(t246.*(t186+t188+t189+t14.*t80.*(t32+t33-t35+t121-t154-t155))+t207.*t229);

dTdK(1) = DTxK;
dTdK(2) = DTyK;