// https://d3js.org/d3-shape/ v1.3.5 Copyright 2019 Mike Bostock
!function(t,n){"object"==typeof exports&&"undefined"!=typeof module?n(exports,require("d3-path")):"function"==typeof define&&define.amd?define(["exports","d3-path"],n):n(t.d3=t.d3||{},t.d3)}(this,function(t,n){"use strict";function i(t){return function(){return t}}var e=Math.abs,s=Math.atan2,h=Math.cos,o=Math.max,_=Math.min,r=Math.sin,a=Math.sqrt,c=1e-12,l=Math.PI,u=l/2,f=2*l;function x(t){return t>=1?u:t<=-1?-u:Math.asin(t)}function y(t){return t.innerRadius}function p(t){return t.outerRadius}function v(t){return t.startAngle}function d(t){return t.endAngle}function T(t){return t&&t.padAngle}function g(t,n,i,e,s,h,_){var r=t-i,c=n-e,l=(_?h:-h)/a(r*r+c*c),u=l*c,f=-l*r,x=t+u,y=n+f,p=i+u,v=e+f,d=(x+p)/2,T=(y+v)/2,g=p-x,b=v-y,w=g*g+b*b,k=s-h,m=x*v-p*y,N=(b<0?-1:1)*a(o(0,k*k*w-m*m)),M=(m*b-g*N)/w,S=(-m*g-b*N)/w,E=(m*b+g*N)/w,A=(-m*g+b*N)/w,P=M-d,C=S-T,O=E-d,R=A-T;return P*P+C*C>O*O+R*R&&(M=E,S=A),{cx:M,cy:S,x01:-u,y01:-f,x11:M*(s/k-1),y11:S*(s/k-1)}}function b(t){this._context=t}function w(t){return new b(t)}function k(t){return t[0]}function m(t){return t[1]}function N(){var t=k,e=m,s=i(!0),h=null,o=w,_=null;function r(i){var r,a,c,l=i.length,u=!1;for(null==h&&(_=o(c=n.path())),r=0;r<=l;++r)!(r<l&&s(a=i[r],r,i))===u&&((u=!u)?_.lineStart():_.lineEnd()),u&&_.point(+t(a,r,i),+e(a,r,i));if(c)return _=null,c+""||null}return r.x=function(n){return arguments.length?(t="function"==typeof n?n:i(+n),r):t},r.y=function(t){return arguments.length?(e="function"==typeof t?t:i(+t),r):e},r.defined=function(t){return arguments.length?(s="function"==typeof t?t:i(!!t),r):s},r.curve=function(t){return arguments.length?(o=t,null!=h&&(_=o(h)),r):o},r.context=function(t){return arguments.length?(null==t?h=_=null:_=o(h=t),r):h},r}function M(){var t=k,e=null,s=i(0),h=m,o=i(!0),_=null,r=w,a=null;function c(i){var c,l,u,f,x,y=i.length,p=!1,v=new Array(y),d=new Array(y);for(null==_&&(a=r(x=n.path())),c=0;c<=y;++c){if(!(c<y&&o(f=i[c],c,i))===p)if(p=!p)l=c,a.areaStart(),a.lineStart();else{for(a.lineEnd(),a.lineStart(),u=c-1;u>=l;--u)a.point(v[u],d[u]);a.lineEnd(),a.areaEnd()}p&&(v[c]=+t(f,c,i),d[c]=+s(f,c,i),a.point(e?+e(f,c,i):v[c],h?+h(f,c,i):d[c]))}if(x)return a=null,x+""||null}function l(){return N().defined(o).curve(r).context(_)}return c.x=function(n){return arguments.length?(t="function"==typeof n?n:i(+n),e=null,c):t},c.x0=function(n){return arguments.length?(t="function"==typeof n?n:i(+n),c):t},c.x1=function(t){return arguments.length?(e=null==t?null:"function"==typeof t?t:i(+t),c):e},c.y=function(t){return arguments.length?(s="function"==typeof t?t:i(+t),h=null,c):s},c.y0=function(t){return arguments.length?(s="function"==typeof t?t:i(+t),c):s},c.y1=function(t){return arguments.length?(h=null==t?null:"function"==typeof t?t:i(+t),c):h},c.lineX0=c.lineY0=function(){return l().x(t).y(s)},c.lineY1=function(){return l().x(t).y(h)},c.lineX1=function(){return l().x(e).y(s)},c.defined=function(t){return arguments.length?(o="function"==typeof t?t:i(!!t),c):o},c.curve=function(t){return arguments.length?(r=t,null!=_&&(a=r(_)),c):r},c.context=function(t){return arguments.length?(null==t?_=a=null:a=r(_=t),c):_},c}function S(t,n){return n<t?-1:n>t?1:n>=t?0:NaN}function E(t){return t}b.prototype={areaStart:function(){this._line=0},areaEnd:function(){this._line=NaN},lineStart:function(){this._point=0},lineEnd:function(){(this._line||0!==this._line&&1===this._point)&&this._context.closePath(),this._line=1-this._line},point:function(t,n){switch(t=+t,n=+n,this._point){case 0:this._point=1,this._line?this._context.lineTo(t,n):this._context.moveTo(t,n);break;case 1:this._point=2;default:this._context.lineTo(t,n)}}};var A=C(w);function P(t){this._curve=t}function C(t){function n(n){return new P(t(n))}return n._curve=t,n}function O(t){var n=t.curve;return t.angle=t.x,delete t.x,t.radius=t.y,delete t.y,t.curve=function(t){return arguments.length?n(C(t)):n()._curve},t}function R(){return O(N().curve(A))}function q(){var t=M().curve(A),n=t.curve,i=t.lineX0,e=t.lineX1,s=t.lineY0,h=t.lineY1;return t.angle=t.x,delete t.x,t.startAngle=t.x0,delete t.x0,t.endAngle=t.x1,delete t.x1,t.radius=t.y,delete t.y,t.innerRadius=t.y0,delete t.y0,t.outerRadius=t.y1,delete t.y1,t.lineStartAngle=function(){return O(i())},delete t.lineX0,t.lineEndAngle=function(){return O(e())},delete t.lineX1,t.lineInnerRadius=function(){return O(s())},delete t.lineY0,t.lineOuterRadius=function(){return O(h())},delete t.lineY1,t.curve=function(t){return arguments.length?n(C(t)):n()._curve},t}function z(t,n){return[(n=+n)*Math.cos(t-=Math.PI/2),n*Math.sin(t)]}P.prototype={areaStart:function(){this._curve.areaStart()},areaEnd:function(){this._curve.areaEnd()},lineStart:function(){this._curve.lineStart()},lineEnd:function(){this._curve.lineEnd()},point:function(t,n){this._curve.point(n*Math.sin(t),n*-Math.cos(t))}};var X=Array.prototype.slice;function Y(t){return t.source}function B(t){return t.target}function I(t){var e=Y,s=B,h=k,o=m,_=null;function r(){var i,r=X.call(arguments),a=e.apply(this,r),c=s.apply(this,r);if(_||(_=i=n.path()),t(_,+h.apply(this,(r[0]=a,r)),+o.apply(this,r),+h.apply(this,(r[0]=c,r)),+o.apply(this,r)),i)return _=null,i+""||null}return r.source=function(t){return arguments.length?(e=t,r):e},r.target=function(t){return arguments.length?(s=t,r):s},r.x=function(t){return arguments.length?(h="function"==typeof t?t:i(+t),r):h},r.y=function(t){return arguments.length?(o="function"==typeof t?t:i(+t),r):o},r.context=function(t){return arguments.length?(_=null==t?null:t,r):_},r}function j(t,n,i,e,s){t.moveTo(n,i),t.bezierCurveTo(n=(n+e)/2,i,n,s,e,s)}function D(t,n,i,e,s){t.moveTo(n,i),t.bezierCurveTo(n,i=(i+s)/2,e,i,e,s)}function L(t,n,i,e,s){var h=z(n,i),o=z(n,i=(i+s)/2),_=z(e,i),r=z(e,s);t.moveTo(h[0],h[1]),t.bezierCurveTo(o[0],o[1],_[0],_[1],r[0],r[1])}var V={draw:function(t,n){var i=Math.sqrt(n/l);t.moveTo(i,0),t.arc(0,0,i,0,f)}},W={draw:function(t,n){var i=Math.sqrt(n/5)/2;t.moveTo(-3*i,-i),t.lineTo(-i,-i),t.lineTo(-i,-3*i),t.lineTo(i,-3*i),t.lineTo(i,-i),t.lineTo(3*i,-i),t.lineTo(3*i,i),t.lineTo(i,i),t.lineTo(i,3*i),t.lineTo(-i,3*i),t.lineTo(-i,i),t.lineTo(-3*i,i),t.closePath()}},H=Math.sqrt(1/3),F=2*H,G={draw:function(t,n){var i=Math.sqrt(n/F),e=i*H;t.moveTo(0,-i),t.lineTo(e,0),t.lineTo(0,i),t.lineTo(-e,0),t.closePath()}},J=Math.sin(l/10)/Math.sin(7*l/10),K=Math.sin(f/10)*J,Q=-Math.cos(f/10)*J,U={draw:function(t,n){var i=Math.sqrt(.8908130915292852*n),e=K*i,s=Q*i;t.moveTo(0,-i),t.lineTo(e,s);for(var h=1;h<5;++h){var o=f*h/5,_=Math.cos(o),r=Math.sin(o);t.lineTo(r*i,-_*i),t.lineTo(_*e-r*s,r*e+_*s)}t.closePath()}},Z={draw:function(t,n){var i=Math.sqrt(n),e=-i/2;t.rect(e,e,i,i)}},$=Math.sqrt(3),tt={draw:function(t,n){var i=-Math.sqrt(n/(3*$));t.moveTo(0,2*i),t.lineTo(-$*i,-i),t.lineTo($*i,-i),t.closePath()}},nt=-.5,it=Math.sqrt(3)/2,et=1/Math.sqrt(12),st=3*(et/2+1),ht={draw:function(t,n){var i=Math.sqrt(n/st),e=i/2,s=i*et,h=e,o=i*et+i,_=-h,r=o;t.moveTo(e,s),t.lineTo(h,o),t.lineTo(_,r),t.lineTo(nt*e-it*s,it*e+nt*s),t.lineTo(nt*h-it*o,it*h+nt*o),t.lineTo(nt*_-it*r,it*_+nt*r),t.lineTo(nt*e+it*s,nt*s-it*e),t.lineTo(nt*h+it*o,nt*o-it*h),t.lineTo(nt*_+it*r,nt*r-it*_),t.closePath()}},ot=[V,W,G,Z,U,tt,ht];function _t(){}function rt(t,n,i){t._context.bezierCurveTo((2*t._x0+t._x1)/3,(2*t._y0+t._y1)/3,(t._x0+2*t._x1)/3,(t._y0+2*t._y1)/3,(t._x0+4*t._x1+n)/6,(t._y0+4*t._y1+i)/6)}function at(t){this._context=t}function ct(t){this._context=t}function lt(t){this._context=t}function ut(t,n){this._basis=new at(t),this._beta=n}at.prototype={areaStart:function(){this._line=0},areaEnd:function(){this._line=NaN},lineStart:function(){this._x0=this._x1=this._y0=this._y1=NaN,this._point=0},lineEnd:function(){switch(this._point){case 3:rt(this,this._x1,this._y1);case 2:this._context.lineTo(this._x1,this._y1)}(this._line||0!==this._line&&1===this._point)&&this._context.closePath(),this._line=1-this._line},point:function(t,n){switch(t=+t,n=+n,this._point){case 0:this._point=1,this._line?this._context.lineTo(t,n):this._context.moveTo(t,n);break;case 1:this._point=2;break;case 2:this._point=3,this._context.lineTo((5*this._x0+this._x1)/6,(5*this._y0+this._y1)/6);default:rt(this,t,n)}this._x0=this._x1,this._x1=t,this._y0=this._y1,this._y1=n}},ct.prototype={areaStart:_t,areaEnd:_t,lineStart:function(){this._x0=this._x1=this._x2=this._x3=this._x4=this._y0=this._y1=this._y2=this._y3=this._y4=NaN,this._point=0},lineEnd:function(){switch(this._point){case 1:this._context.moveTo(this._x2,this._y2),this._context.closePath();break;case 2:this._context.moveTo((this._x2+2*this._x3)/3,(this._y2+2*this._y3)/3),this._context.lineTo((this._x3+2*this._x2)/3,(this._y3+2*this._y2)/3),this._context.closePath();break;case 3:this.point(this._x2,this._y2),this.point(this._x3,this._y3),this.point(this._x4,this._y4)}},point:function(t,n){switch(t=+t,n=+n,this._point){case 0:this._point=1,this._x2=t,this._y2=n;break;case 1:this._point=2,this._x3=t,this._y3=n;break;case 2:this._point=3,this._x4=t,this._y4=n,this._context.moveTo((this._x0+4*this._x1+t)/6,(this._y0+4*this._y1+n)/6);break;default:rt(this,t,n)}this._x0=this._x1,this._x1=t,this._y0=this._y1,this._y1=n}},lt.prototype={areaStart:function(){this._line=0},areaEnd:function(){this._line=NaN},lineStart:function(){this._x0=this._x1=this._y0=this._y1=NaN,this._point=0},lineEnd:function(){(this._line||0!==this._line&&3===this._point)&&this._context.closePath(),this._line=1-this._line},point:function(t,n){switch(t=+t,n=+n,this._point){case 0:this._point=1;break;case 1:this._point=2;break;case 2:this._point=3;var i=(this._x0+4*this._x1+t)/6,e=(this._y0+4*this._y1+n)/6;this._line?this._context.lineTo(i,e):this._context.moveTo(i,e);break;case 3:this._point=4;default:rt(this,t,n)}this._x0=this._x1,this._x1=t,this._y0=this._y1,this._y1=n}},ut.prototype={lineStart:function(){this._x=[],this._y=[],this._basis.lineStart()},lineEnd:function(){var t=this._x,n=this._y,i=t.length-1;if(i>0)for(var e,s=t[0],h=n[0],o=t[i]-s,_=n[i]-h,r=-1;++r<=i;)e=r/i,this._basis.point(this._beta*t[r]+(1-this._beta)*(s+e*o),this._beta*n[r]+(1-this._beta)*(h+e*_));this._x=this._y=null,this._basis.lineEnd()},point:function(t,n){this._x.push(+t),this._y.push(+n)}};var ft=function t(n){function i(t){return 1===n?new at(t):new ut(t,n)}return i.beta=function(n){return t(+n)},i}(.85);function xt(t,n,i){t._context.bezierCurveTo(t._x1+t._k*(t._x2-t._x0),t._y1+t._k*(t._y2-t._y0),t._x2+t._k*(t._x1-n),t._y2+t._k*(t._y1-i),t._x2,t._y2)}function yt(t,n){this._context=t,this._k=(1-n)/6}yt.prototype={areaStart:function(){this._line=0},areaEnd:function(){this._line=NaN},lineStart:function(){this._x0=this._x1=this._x2=this._y0=this._y1=this._y2=NaN,this._point=0},lineEnd:function(){switch(this._point){case 2:this._context.lineTo(this._x2,this._y2);break;case 3:xt(this,this._x1,this._y1)}(this._line||0!==this._line&&1===this._point)&&this._context.closePath(),this._line=1-this._line},point:function(t,n){switch(t=+t,n=+n,this._point){case 0:this._point=1,this._line?this._context.lineTo(t,n):this._context.moveTo(t,n);break;case 1:this._point=2,this._x1=t,this._y1=n;break;case 2:this._point=3;default:xt(this,t,n)}this._x0=this._x1,this._x1=this._x2,this._x2=t,this._y0=this._y1,this._y1=this._y2,this._y2=n}};var pt=function t(n){function i(t){return new yt(t,n)}return i.tension=function(n){return t(+n)},i}(0);function vt(t,n){this._context=t,this._k=(1-n)/6}vt.prototype={areaStart:_t,areaEnd:_t,lineStart:function(){this._x0=this._x1=this._x2=this._x3=this._x4=this._x5=this._y0=this._y1=this._y2=this._y3=this._y4=this._y5=NaN,this._point=0},lineEnd:function(){switch(this._point){case 1:this._context.moveTo(this._x3,this._y3),this._context.closePath();break;case 2:this._context.lineTo(this._x3,this._y3),this._context.closePath();break;case 3:this.point(this._x3,this._y3),this.point(this._x4,this._y4),this.point(this._x5,this._y5)}},point:function(t,n){switch(t=+t,n=+n,this._point){case 0:this._point=1,this._x3=t,this._y3=n;break;case 1:this._point=2,this._context.moveTo(this._x4=t,this._y4=n);break;case 2:this._point=3,this._x5=t,this._y5=n;break;default:xt(this,t,n)}this._x0=this._x1,this._x1=this._x2,this._x2=t,this._y0=this._y1,this._y1=this._y2,this._y2=n}};var dt=function t(n){function i(t){return new vt(t,n)}return i.tension=function(n){return t(+n)},i}(0);function Tt(t,n){this._context=t,this._k=(1-n)/6}Tt.prototype={areaStart:function(){this._line=0},areaEnd:function(){this._line=NaN},lineStart:function(){this._x0=this._x1=this._x2=this._y0=this._y1=this._y2=NaN,this._point=0},lineEnd:function(){(this._line||0!==this._line&&3===this._point)&&this._context.closePath(),this._line=1-this._line},point:function(t,n){switch(t=+t,n=+n,this._point){case 0:this._point=1;break;case 1:this._point=2;break;case 2:this._point=3,this._line?this._context.lineTo(this._x2,this._y2):this._context.moveTo(this._x2,this._y2);break;case 3:this._point=4;default:xt(this,t,n)}this._x0=this._x1,this._x1=this._x2,this._x2=t,this._y0=this._y1,this._y1=this._y2,this._y2=n}};var gt=function t(n){function i(t){return new Tt(t,n)}return i.tension=function(n){return t(+n)},i}(0);function bt(t,n,i){var e=t._x1,s=t._y1,h=t._x2,o=t._y2;if(t._l01_a>c){var _=2*t._l01_2a+3*t._l01_a*t._l12_a+t._l12_2a,r=3*t._l01_a*(t._l01_a+t._l12_a);e=(e*_-t._x0*t._l12_2a+t._x2*t._l01_2a)/r,s=(s*_-t._y0*t._l12_2a+t._y2*t._l01_2a)/r}if(t._l23_a>c){var a=2*t._l23_2a+3*t._l23_a*t._l12_a+t._l12_2a,l=3*t._l23_a*(t._l23_a+t._l12_a);h=(h*a+t._x1*t._l23_2a-n*t._l12_2a)/l,o=(o*a+t._y1*t._l23_2a-i*t._l12_2a)/l}t._context.bezierCurveTo(e,s,h,o,t._x2,t._y2)}function wt(t,n){this._context=t,this._alpha=n}wt.prototype={areaStart:function(){this._line=0},areaEnd:function(){this._line=NaN},lineStart:function(){this._x0=this._x1=this._x2=this._y0=this._y1=this._y2=NaN,this._l01_a=this._l12_a=this._l23_a=this._l01_2a=this._l12_2a=this._l23_2a=this._point=0},lineEnd:function(){switch(this._point){case 2:this._context.lineTo(this._x2,this._y2);break;case 3:this.point(this._x2,this._y2)}(this._line||0!==this._line&&1===this._point)&&this._context.closePath(),this._line=1-this._line},point:function(t,n){if(t=+t,n=+n,this._point){var i=this._x2-t,e=this._y2-n;this._l23_a=Math.sqrt(this._l23_2a=Math.pow(i*i+e*e,this._alpha))}switch(this._point){case 0:this._point=1,this._line?this._context.lineTo(t,n):this._context.moveTo(t,n);break;case 1:this._point=2;break;case 2:this._point=3;default:bt(this,t,n)}this._l01_a=this._l12_a,this._l12_a=this._l23_a,this._l01_2a=this._l12_2a,this._l12_2a=this._l23_2a,this._x0=this._x1,this._x1=this._x2,this._x2=t,this._y0=this._y1,this._y1=this._y2,this._y2=n}};var kt=function t(n){function i(t){return n?new wt(t,n):new yt(t,0)}return i.alpha=function(n){return t(+n)},i}(.5);function mt(t,n){this._context=t,this._alpha=n}mt.prototype={areaStart:_t,areaEnd:_t,lineStart:function(){this._x0=this._x1=this._x2=this._x3=this._x4=this._x5=this._y0=this._y1=this._y2=this._y3=this._y4=this._y5=NaN,this._l01_a=this._l12_a=this._l23_a=this._l01_2a=this._l12_2a=this._l23_2a=this._point=0},lineEnd:function(){switch(this._point){case 1:this._context.moveTo(this._x3,this._y3),this._context.closePath();break;case 2:this._context.lineTo(this._x3,this._y3),this._context.closePath();break;case 3:this.point(this._x3,this._y3),this.point(this._x4,this._y4),this.point(this._x5,this._y5)}},point:function(t,n){if(t=+t,n=+n,this._point){var i=this._x2-t,e=this._y2-n;this._l23_a=Math.sqrt(this._l23_2a=Math.pow(i*i+e*e,this._alpha))}switch(this._point){case 0:this._point=1,this._x3=t,this._y3=n;break;case 1:this._point=2,this._context.moveTo(this._x4=t,this._y4=n);break;case 2:this._point=3,this._x5=t,this._y5=n;break;default:bt(this,t,n)}this._l01_a=this._l12_a,this._l12_a=this._l23_a,this._l01_2a=this._l12_2a,this._l12_2a=this._l23_2a,this._x0=this._x1,this._x1=this._x2,this._x2=t,this._y0=this._y1,this._y1=this._y2,this._y2=n}};var Nt=function t(n){function i(t){return n?new mt(t,n):new vt(t,0)}return i.alpha=function(n){return t(+n)},i}(.5);function Mt(t,n){this._context=t,this._alpha=n}Mt.prototype={areaStart:function(){this._line=0},areaEnd:function(){this._line=NaN},lineStart:function(){this._x0=this._x1=this._x2=this._y0=this._y1=this._y2=NaN,this._l01_a=this._l12_a=this._l23_a=this._l01_2a=this._l12_2a=this._l23_2a=this._point=0},lineEnd:function(){(this._line||0!==this._line&&3===this._point)&&this._context.closePath(),this._line=1-this._line},point:function(t,n){if(t=+t,n=+n,this._point){var i=this._x2-t,e=this._y2-n;this._l23_a=Math.sqrt(this._l23_2a=Math.pow(i*i+e*e,this._alpha))}switch(this._point){case 0:this._point=1;break;case 1:this._point=2;break;case 2:this._point=3,this._line?this._context.lineTo(this._x2,this._y2):this._context.moveTo(this._x2,this._y2);break;case 3:this._point=4;default:bt(this,t,n)}this._l01_a=this._l12_a,this._l12_a=this._l23_a,this._l01_2a=this._l12_2a,this._l12_2a=this._l23_2a,this._x0=this._x1,this._x1=this._x2,this._x2=t,this._y0=this._y1,this._y1=this._y2,this._y2=n}};var St=function t(n){function i(t){return n?new Mt(t,n):new Tt(t,0)}return i.alpha=function(n){return t(+n)},i}(.5);function Et(t){this._context=t}function At(t){return t<0?-1:1}function Pt(t,n,i){var e=t._x1-t._x0,s=n-t._x1,h=(t._y1-t._y0)/(e||s<0&&-0),o=(i-t._y1)/(s||e<0&&-0),_=(h*s+o*e)/(e+s);return(At(h)+At(o))*Math.min(Math.abs(h),Math.abs(o),.5*Math.abs(_))||0}function Ct(t,n){var i=t._x1-t._x0;return i?(3*(t._y1-t._y0)/i-n)/2:n}function Ot(t,n,i){var e=t._x0,s=t._y0,h=t._x1,o=t._y1,_=(h-e)/3;t._context.bezierCurveTo(e+_,s+_*n,h-_,o-_*i,h,o)}function Rt(t){this._context=t}function qt(t){this._context=new zt(t)}function zt(t){this._context=t}function Xt(t){this._context=t}function Yt(t){var n,i,e=t.length-1,s=new Array(e),h=new Array(e),o=new Array(e);for(s[0]=0,h[0]=2,o[0]=t[0]+2*t[1],n=1;n<e-1;++n)s[n]=1,h[n]=4,o[n]=4*t[n]+2*t[n+1];for(s[e-1]=2,h[e-1]=7,o[e-1]=8*t[e-1]+t[e],n=1;n<e;++n)i=s[n]/h[n-1],h[n]-=i,o[n]-=i*o[n-1];for(s[e-1]=o[e-1]/h[e-1],n=e-2;n>=0;--n)s[n]=(o[n]-s[n+1])/h[n];for(h[e-1]=(t[e]+s[e-1])/2,n=0;n<e-1;++n)h[n]=2*t[n+1]-s[n+1];return[s,h]}function Bt(t,n){this._context=t,this._t=n}function It(t,n){if((s=t.length)>1)for(var i,e,s,h=1,o=t[n[0]],_=o.length;h<s;++h)for(e=o,o=t[n[h]],i=0;i<_;++i)o[i][1]+=o[i][0]=isNaN(e[i][1])?e[i][0]:e[i][1]}function jt(t){for(var n=t.length,i=new Array(n);--n>=0;)i[n]=n;return i}function Dt(t,n){return t[n]}function Lt(t){var n=t.map(Vt);return jt(t).sort(function(t,i){return n[t]-n[i]})}function Vt(t){for(var n,i=-1,e=0,s=t.length,h=-1/0;++i<s;)(n=+t[i][1])>h&&(h=n,e=i);return e}function Wt(t){var n=t.map(Ht);return jt(t).sort(function(t,i){return n[t]-n[i]})}function Ht(t){for(var n,i=0,e=-1,s=t.length;++e<s;)(n=+t[e][1])&&(i+=n);return i}Et.prototype={areaStart:_t,areaEnd:_t,lineStart:function(){this._point=0},lineEnd:function(){this._point&&this._context.closePath()},point:function(t,n){t=+t,n=+n,this._point?this._context.lineTo(t,n):(this._point=1,this._context.moveTo(t,n))}},Rt.prototype={areaStart:function(){this._line=0},areaEnd:function(){this._line=NaN},lineStart:function(){this._x0=this._x1=this._y0=this._y1=this._t0=NaN,this._point=0},lineEnd:function(){switch(this._point){case 2:this._context.lineTo(this._x1,this._y1);break;case 3:Ot(this,this._t0,Ct(this,this._t0))}(this._line||0!==this._line&&1===this._point)&&this._context.closePath(),this._line=1-this._line},point:function(t,n){var i=NaN;if(n=+n,(t=+t)!==this._x1||n!==this._y1){switch(this._point){case 0:this._point=1,this._line?this._context.lineTo(t,n):this._context.moveTo(t,n);break;case 1:this._point=2;break;case 2:this._point=3,Ot(this,Ct(this,i=Pt(this,t,n)),i);break;default:Ot(this,this._t0,i=Pt(this,t,n))}this._x0=this._x1,this._x1=t,this._y0=this._y1,this._y1=n,this._t0=i}}},(qt.prototype=Object.create(Rt.prototype)).point=function(t,n){Rt.prototype.point.call(this,n,t)},zt.prototype={moveTo:function(t,n){this._context.moveTo(n,t)},closePath:function(){this._context.closePath()},lineTo:function(t,n){this._context.lineTo(n,t)},bezierCurveTo:function(t,n,i,e,s,h){this._context.bezierCurveTo(n,t,e,i,h,s)}},Xt.prototype={areaStart:function(){this._line=0},areaEnd:function(){this._line=NaN},lineStart:function(){this._x=[],this._y=[]},lineEnd:function(){var t=this._x,n=this._y,i=t.length;if(i)if(this._line?this._context.lineTo(t[0],n[0]):this._context.moveTo(t[0],n[0]),2===i)this._context.lineTo(t[1],n[1]);else for(var e=Yt(t),s=Yt(n),h=0,o=1;o<i;++h,++o)this._context.bezierCurveTo(e[0][h],s[0][h],e[1][h],s[1][h],t[o],n[o]);(this._line||0!==this._line&&1===i)&&this._context.closePath(),this._line=1-this._line,this._x=this._y=null},point:function(t,n){this._x.push(+t),this._y.push(+n)}},Bt.prototype={areaStart:function(){this._line=0},areaEnd:function(){this._line=NaN},lineStart:function(){this._x=this._y=NaN,this._point=0},lineEnd:function(){0<this._t&&this._t<1&&2===this._point&&this._context.lineTo(this._x,this._y),(this._line||0!==this._line&&1===this._point)&&this._context.closePath(),this._line>=0&&(this._t=1-this._t,this._line=1-this._line)},point:function(t,n){switch(t=+t,n=+n,this._point){case 0:this._point=1,this._line?this._context.lineTo(t,n):this._context.moveTo(t,n);break;case 1:this._point=2;default:if(this._t<=0)this._context.lineTo(this._x,n),this._context.lineTo(t,n);else{var i=this._x*(1-this._t)+t*this._t;this._context.lineTo(i,this._y),this._context.lineTo(i,n)}}this._x=t,this._y=n}},t.arc=function(){var t=y,o=p,b=i(0),w=null,k=v,m=d,N=T,M=null;function S(){var i,y,p,v=+t.apply(this,arguments),d=+o.apply(this,arguments),T=k.apply(this,arguments)-u,S=m.apply(this,arguments)-u,E=e(S-T),A=S>T;if(M||(M=i=n.path()),d<v&&(y=d,d=v,v=y),d>c)if(E>f-c)M.moveTo(d*h(T),d*r(T)),M.arc(0,0,d,T,S,!A),v>c&&(M.moveTo(v*h(S),v*r(S)),M.arc(0,0,v,S,T,A));else{var P,C,O=T,R=S,q=T,z=S,X=E,Y=E,B=N.apply(this,arguments)/2,I=B>c&&(w?+w.apply(this,arguments):a(v*v+d*d)),j=_(e(d-v)/2,+b.apply(this,arguments)),D=j,L=j;if(I>c){var V=x(I/v*r(B)),W=x(I/d*r(B));(X-=2*V)>c?(q+=V*=A?1:-1,z-=V):(X=0,q=z=(T+S)/2),(Y-=2*W)>c?(O+=W*=A?1:-1,R-=W):(Y=0,O=R=(T+S)/2)}var H=d*h(O),F=d*r(O),G=v*h(z),J=v*r(z);if(j>c){var K,Q=d*h(R),U=d*r(R),Z=v*h(q),$=v*r(q);if(E<l&&(K=function(t,n,i,e,s,h,o,_){var r=i-t,a=e-n,l=o-s,u=_-h,f=u*r-l*a;if(!(f*f<c))return[t+(f=(l*(n-h)-u*(t-s))/f)*r,n+f*a]}(H,F,Z,$,Q,U,G,J))){var tt=H-K[0],nt=F-K[1],it=Q-K[0],et=U-K[1],st=1/r(((p=(tt*it+nt*et)/(a(tt*tt+nt*nt)*a(it*it+et*et)))>1?0:p<-1?l:Math.acos(p))/2),ht=a(K[0]*K[0]+K[1]*K[1]);D=_(j,(v-ht)/(st-1)),L=_(j,(d-ht)/(st+1))}}Y>c?L>c?(P=g(Z,$,H,F,d,L,A),C=g(Q,U,G,J,d,L,A),M.moveTo(P.cx+P.x01,P.cy+P.y01),L<j?M.arc(P.cx,P.cy,L,s(P.y01,P.x01),s(C.y01,C.x01),!A):(M.arc(P.cx,P.cy,L,s(P.y01,P.x01),s(P.y11,P.x11),!A),M.arc(0,0,d,s(P.cy+P.y11,P.cx+P.x11),s(C.cy+C.y11,C.cx+C.x11),!A),M.arc(C.cx,C.cy,L,s(C.y11,C.x11),s(C.y01,C.x01),!A))):(M.moveTo(H,F),M.arc(0,0,d,O,R,!A)):M.moveTo(H,F),v>c&&X>c?D>c?(P=g(G,J,Q,U,v,-D,A),C=g(H,F,Z,$,v,-D,A),M.lineTo(P.cx+P.x01,P.cy+P.y01),D<j?M.arc(P.cx,P.cy,D,s(P.y01,P.x01),s(C.y01,C.x01),!A):(M.arc(P.cx,P.cy,D,s(P.y01,P.x01),s(P.y11,P.x11),!A),M.arc(0,0,v,s(P.cy+P.y11,P.cx+P.x11),s(C.cy+C.y11,C.cx+C.x11),A),M.arc(C.cx,C.cy,D,s(C.y11,C.x11),s(C.y01,C.x01),!A))):M.arc(0,0,v,z,q,A):M.lineTo(G,J)}else M.moveTo(0,0);if(M.closePath(),i)return M=null,i+""||null}return S.centroid=function(){var n=(+t.apply(this,arguments)+ +o.apply(this,arguments))/2,i=(+k.apply(this,arguments)+ +m.apply(this,arguments))/2-l/2;return[h(i)*n,r(i)*n]},S.innerRadius=function(n){return arguments.length?(t="function"==typeof n?n:i(+n),S):t},S.outerRadius=function(t){return arguments.length?(o="function"==typeof t?t:i(+t),S):o},S.cornerRadius=function(t){return arguments.length?(b="function"==typeof t?t:i(+t),S):b},S.padRadius=function(t){return arguments.length?(w=null==t?null:"function"==typeof t?t:i(+t),S):w},S.startAngle=function(t){return arguments.length?(k="function"==typeof t?t:i(+t),S):k},S.endAngle=function(t){return arguments.length?(m="function"==typeof t?t:i(+t),S):m},S.padAngle=function(t){return arguments.length?(N="function"==typeof t?t:i(+t),S):N},S.context=function(t){return arguments.length?(M=null==t?null:t,S):M},S},t.area=M,t.line=N,t.pie=function(){var t=E,n=S,e=null,s=i(0),h=i(f),o=i(0);function _(i){var _,r,a,c,l,u=i.length,x=0,y=new Array(u),p=new Array(u),v=+s.apply(this,arguments),d=Math.min(f,Math.max(-f,h.apply(this,arguments)-v)),T=Math.min(Math.abs(d)/u,o.apply(this,arguments)),g=T*(d<0?-1:1);for(_=0;_<u;++_)(l=p[y[_]=_]=+t(i[_],_,i))>0&&(x+=l);for(null!=n?y.sort(function(t,i){return n(p[t],p[i])}):null!=e&&y.sort(function(t,n){return e(i[t],i[n])}),_=0,a=x?(d-u*g)/x:0;_<u;++_,v=c)r=y[_],c=v+((l=p[r])>0?l*a:0)+g,p[r]={data:i[r],index:_,value:l,startAngle:v,endAngle:c,padAngle:T};return p}return _.value=function(n){return arguments.length?(t="function"==typeof n?n:i(+n),_):t},_.sortValues=function(t){return arguments.length?(n=t,e=null,_):n},_.sort=function(t){return arguments.length?(e=t,n=null,_):e},_.startAngle=function(t){return arguments.length?(s="function"==typeof t?t:i(+t),_):s},_.endAngle=function(t){return arguments.length?(h="function"==typeof t?t:i(+t),_):h},_.padAngle=function(t){return arguments.length?(o="function"==typeof t?t:i(+t),_):o},_},t.areaRadial=q,t.radialArea=q,t.lineRadial=R,t.radialLine=R,t.pointRadial=z,t.linkHorizontal=function(){return I(j)},t.linkVertical=function(){return I(D)},t.linkRadial=function(){var t=I(L);return t.angle=t.x,delete t.x,t.radius=t.y,delete t.y,t},t.symbol=function(){var t=i(V),e=i(64),s=null;function h(){var i;if(s||(s=i=n.path()),t.apply(this,arguments).draw(s,+e.apply(this,arguments)),i)return s=null,i+""||null}return h.type=function(n){return arguments.length?(t="function"==typeof n?n:i(n),h):t},h.size=function(t){return arguments.length?(e="function"==typeof t?t:i(+t),h):e},h.context=function(t){return arguments.length?(s=null==t?null:t,h):s},h},t.symbols=ot,t.symbolCircle=V,t.symbolCross=W,t.symbolDiamond=G,t.symbolSquare=Z,t.symbolStar=U,t.symbolTriangle=tt,t.symbolWye=ht,t.curveBasisClosed=function(t){return new ct(t)},t.curveBasisOpen=function(t){return new lt(t)},t.curveBasis=function(t){return new at(t)},t.curveBundle=ft,t.curveCardinalClosed=dt,t.curveCardinalOpen=gt,t.curveCardinal=pt,t.curveCatmullRomClosed=Nt,t.curveCatmullRomOpen=St,t.curveCatmullRom=kt,t.curveLinearClosed=function(t){return new Et(t)},t.curveLinear=w,t.curveMonotoneX=function(t){return new Rt(t)},t.curveMonotoneY=function(t){return new qt(t)},t.curveNatural=function(t){return new Xt(t)},t.curveStep=function(t){return new Bt(t,.5)},t.curveStepAfter=function(t){return new Bt(t,1)},t.curveStepBefore=function(t){return new Bt(t,0)},t.stack=function(){var t=i([]),n=jt,e=It,s=Dt;function h(i){var h,o,_=t.apply(this,arguments),r=i.length,a=_.length,c=new Array(a);for(h=0;h<a;++h){for(var l,u=_[h],f=c[h]=new Array(r),x=0;x<r;++x)f[x]=l=[0,+s(i[x],u,x,i)],l.data=i[x];f.key=u}for(h=0,o=n(c);h<a;++h)c[o[h]].index=h;return e(c,o),c}return h.keys=function(n){return arguments.length?(t="function"==typeof n?n:i(X.call(n)),h):t},h.value=function(t){return arguments.length?(s="function"==typeof t?t:i(+t),h):s},h.order=function(t){return arguments.length?(n=null==t?jt:"function"==typeof t?t:i(X.call(t)),h):n},h.offset=function(t){return arguments.length?(e=null==t?It:t,h):e},h},t.stackOffsetExpand=function(t,n){if((e=t.length)>0){for(var i,e,s,h=0,o=t[0].length;h<o;++h){for(s=i=0;i<e;++i)s+=t[i][h][1]||0;if(s)for(i=0;i<e;++i)t[i][h][1]/=s}It(t,n)}},t.stackOffsetDiverging=function(t,n){if((_=t.length)>0)for(var i,e,s,h,o,_,r=0,a=t[n[0]].length;r<a;++r)for(h=o=0,i=0;i<_;++i)(s=(e=t[n[i]][r])[1]-e[0])>=0?(e[0]=h,e[1]=h+=s):s<0?(e[1]=o,e[0]=o+=s):e[0]=h},t.stackOffsetNone=It,t.stackOffsetSilhouette=function(t,n){if((i=t.length)>0){for(var i,e=0,s=t[n[0]],h=s.length;e<h;++e){for(var o=0,_=0;o<i;++o)_+=t[o][e][1]||0;s[e][1]+=s[e][0]=-_/2}It(t,n)}},t.stackOffsetWiggle=function(t,n){if((s=t.length)>0&&(e=(i=t[n[0]]).length)>0){for(var i,e,s,h=0,o=1;o<e;++o){for(var _=0,r=0,a=0;_<s;++_){for(var c=t[n[_]],l=c[o][1]||0,u=(l-(c[o-1][1]||0))/2,f=0;f<_;++f){var x=t[n[f]];u+=(x[o][1]||0)-(x[o-1][1]||0)}r+=l,a+=u*l}i[o-1][1]+=i[o-1][0]=h,r&&(h-=a/r)}i[o-1][1]+=i[o-1][0]=h,It(t,n)}},t.stackOrderAppearance=Lt,t.stackOrderAscending=Wt,t.stackOrderDescending=function(t){return Wt(t).reverse()},t.stackOrderInsideOut=function(t){var n,i,e=t.length,s=t.map(Ht),h=Lt(t),o=0,_=0,r=[],a=[];for(n=0;n<e;++n)i=h[n],o<_?(o+=s[i],r.push(i)):(_+=s[i],a.push(i));return a.reverse().concat(r)},t.stackOrderNone=jt,t.stackOrderReverse=function(t){return jt(t).reverse()},Object.defineProperty(t,"__esModule",{value:!0})});
