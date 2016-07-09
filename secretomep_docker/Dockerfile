FROM ubuntu:14.04
MAINTAINER Andrew Lonsdale, andrew.lonsdale@lonsbio.com.au

RUN apt-get update -y
RUN apt-get install -y tcsh gnuplot netpbm gawk

ADD secretomep-1.0 /secretomep-1.0/
ADD sub_programs /sub_programs/


ADD run /

ENTRYPOINT ["/run"]





