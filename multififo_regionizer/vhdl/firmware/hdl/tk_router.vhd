library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.regionizer_data.all;

entity tk_router is
    port(
      ap_clk   : IN STD_LOGIC;
      enabled  : IN STD_LOGIC;
      newevent : IN STD_LOGIC;
      links_in      : IN particles(NTKSECTORS*NTKFIBERS-1 downto 0);
      fifo_in       : OUT particles(NTKSECTORS*NTKFIFOS-1 downto 0);
      fifo_in_write : OUT std_logic_vector(NTKSECTORS*NTKFIFOS-1 downto 0);
      fifo_in_roll  : OUT std_logic_vector(NTKSECTORS*NTKFIFOS-1 downto 0)
    );
end tk_router;

architecture Behavioral of tk_router is
begin

     reg_input_first : entity work.tk_router_element
            port map(ap_clk => ap_clk, 
                     enabled => enabled,
                     newevent => newevent,
                     links_in => links_in(NTKFIBERS-1 downto 0),
                     fifo_same => fifo_in(     0      *NTKFIFOS+1*NTKFIBERS-1 downto      0      *NTKFIFOS),
                     fifo_next => fifo_in(     1      *NTKFIFOS+2*NTKFIBERS-1 downto      1      *NTKFIFOS+1*NTKFIBERS),
                     fifo_prev => fifo_in((NTKSECTORS-1)*NTKFIFOS+3*NTKFIBERS-1 downto (NTKSECTORS-1)*NTKFIFOS+2*NTKFIBERS),
                     fifo_same_write => fifo_in_write(     0      *NTKFIFOS+1*NTKFIBERS-1 downto      0      *NTKFIFOS),
                     fifo_next_write => fifo_in_write(     1      *NTKFIFOS+2*NTKFIBERS-1 downto      1      *NTKFIFOS+1*NTKFIBERS),
                     fifo_prev_write => fifo_in_write((NTKSECTORS-1)*NTKFIFOS+3*NTKFIBERS-1 downto (NTKSECTORS-1)*NTKFIFOS+2*NTKFIBERS),
                     fifo_same_roll  => fifo_in_roll (     0      *NTKFIFOS+1*NTKFIBERS-1 downto      0      *NTKFIFOS),
                     fifo_next_roll  => fifo_in_roll (     1      *NTKFIFOS+2*NTKFIBERS-1 downto      1      *NTKFIFOS+1*NTKFIBERS),
                     fifo_prev_roll  => fifo_in_roll ((NTKSECTORS-1)*NTKFIFOS+3*NTKFIBERS-1 downto (NTKSECTORS-1)*NTKFIFOS+2*NTKFIBERS)
                 );
    gen_inputs: for isec in NTKSECTORS-2 downto 1 generate
         reg_input : entity work.tk_router_element
                port map(ap_clk => ap_clk, 
                         enabled => enabled,
                         newevent => newevent,
                         links_in => links_in((isec+1)*NTKFIBERS-1 downto isec*NTKFIBERS),
                         fifo_same => fifo_in( isec   *NTKFIFOS+1*NTKFIBERS-1 downto  isec   *NTKFIFOS),
                         fifo_next => fifo_in((isec+1)*NTKFIFOS+2*NTKFIBERS-1 downto (isec+1)*NTKFIFOS+1*NTKFIBERS),
                         fifo_prev => fifo_in((isec-1)*NTKFIFOS+3*NTKFIBERS-1 downto (isec-1)*NTKFIFOS+2*NTKFIBERS),
                         fifo_same_write => fifo_in_write( isec   *NTKFIFOS+1*NTKFIBERS-1 downto  isec   *NTKFIFOS),
                         fifo_next_write => fifo_in_write((isec+1)*NTKFIFOS+2*NTKFIBERS-1 downto (isec+1)*NTKFIFOS+1*NTKFIBERS),
                         fifo_prev_write => fifo_in_write((isec-1)*NTKFIFOS+3*NTKFIBERS-1 downto (isec-1)*NTKFIFOS+2*NTKFIBERS),
                         fifo_same_roll  => fifo_in_roll ( isec   *NTKFIFOS+1*NTKFIBERS-1 downto  isec   *NTKFIFOS),
                         fifo_next_roll  => fifo_in_roll ((isec+1)*NTKFIFOS+2*NTKFIBERS-1 downto (isec+1)*NTKFIFOS+1*NTKFIBERS),
                         fifo_prev_roll  => fifo_in_roll ((isec-1)*NTKFIFOS+3*NTKFIBERS-1 downto (isec-1)*NTKFIFOS+2*NTKFIBERS)
                     );
        end generate gen_inputs;

     reg_input_last : entity work.tk_router_element
            port map(ap_clk => ap_clk, 
                     enabled => enabled,
                     newevent => newevent,
                     links_in => links_in(NTKSECTORS*NTKFIBERS-1 downto (NTKSECTORS-1)*NTKFIBERS),
                     fifo_same => fifo_in((NTKSECTORS-1)*NTKFIFOS+1*NTKFIBERS-1 downto (NTKSECTORS-1)*NTKFIFOS),
                     fifo_next => fifo_in(     0      *NTKFIFOS+2*NTKFIBERS-1 downto      0      *NTKFIFOS+1*NTKFIBERS),
                     fifo_prev => fifo_in((NTKSECTORS-2)*NTKFIFOS+3*NTKFIBERS-1 downto (NTKSECTORS-2)*NTKFIFOS+2*NTKFIBERS),
                     fifo_same_write => fifo_in_write((NTKSECTORS-1)*NTKFIFOS+1*NTKFIBERS-1 downto (NTKSECTORS-1)*NTKFIFOS),
                     fifo_next_write => fifo_in_write(     0      *NTKFIFOS+2*NTKFIBERS-1 downto      0      *NTKFIFOS+1*NTKFIBERS),
                     fifo_prev_write => fifo_in_write((NTKSECTORS-2)*NTKFIFOS+3*NTKFIBERS-1 downto (NTKSECTORS-2)*NTKFIFOS+2*NTKFIBERS),
                     fifo_same_roll  => fifo_in_roll ((NTKSECTORS-1)*NTKFIFOS+1*NTKFIBERS-1 downto (NTKSECTORS-1)*NTKFIFOS),
                     fifo_next_roll  => fifo_in_roll (     0      *NTKFIFOS+2*NTKFIBERS-1 downto      0      *NTKFIFOS+1*NTKFIBERS),
                     fifo_prev_roll  => fifo_in_roll ((NTKSECTORS-2)*NTKFIFOS+3*NTKFIBERS-1 downto (NTKSECTORS-2)*NTKFIFOS+2*NTKFIBERS)
                 );

end Behavioral;
