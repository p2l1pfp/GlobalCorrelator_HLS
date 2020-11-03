library std;
use std.textio.all;
use std.env.all;
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.std_logic_textio.all;

use work.regionizer_data.all;


entity testbench is
--  Port ( );
end testbench;

architecture Behavioral of testbench is
    constant NREGIONS : natural := NTKSECTORS;

    signal clk : std_logic := '0';
    signal rst : std_logic := '0';
    signal start, ready, idle, done : std_logic;
    signal newevent, newevent_out : std_logic;

    signal p_in:  w64s(NTKSECTORS*NTKFIBERS-1 downto 0);
    signal p_out: w64s(NREGIONS-1         downto 0);

    file Fi : text open read_mode is "input-tk.txt";
    file Fo : text open write_mode is "output-tk-vhdl_tb.txt";


begin
    clk  <= not clk after 1.25 ns;
    
    uut : entity work.tk_regionizer
        port map(ap_clk => clk, 
                 ap_rst => rst, 
                 ap_start => start,
                 ap_ready => ready,
                 ap_idle =>  idle,
                 ap_done => done,
                 tracks_in_0_0_V => p_in( 0),
                 tracks_in_0_1_V => p_in( 1),
                 tracks_in_1_0_V => p_in( 2),
                 tracks_in_1_1_V => p_in( 3), 
                 tracks_in_2_0_V => p_in( 4),
                 tracks_in_2_1_V => p_in( 5),
                 tracks_in_3_0_V => p_in( 6),
                 tracks_in_3_1_V => p_in( 7),
                 tracks_in_4_0_V => p_in( 8),
                 tracks_in_4_1_V => p_in( 9), 
                 tracks_in_5_0_V => p_in(10),
                 tracks_in_5_1_V => p_in(11),
                 tracks_in_6_0_V => p_in(12),
                 tracks_in_6_1_V => p_in(13),
                 tracks_in_7_0_V => p_in(14),
                 tracks_in_7_1_V => p_in(15), 
                 tracks_in_8_0_V => p_in(16),
                 tracks_in_8_1_V => p_in(17),
                 tracks_out_0_V => p_out(0),
                 tracks_out_1_V => p_out(1),
                 tracks_out_2_V => p_out(2),
                 tracks_out_3_V => p_out(3),
                 tracks_out_4_V => p_out(4),
                 tracks_out_5_V => p_out(5),
                 tracks_out_6_V => p_out(6),
                 tracks_out_7_V => p_out(7),
                 tracks_out_8_V => p_out(8),
                 newevent => newevent,
                 newevent_out => newevent_out
             );
   

    runit : process 
        variable remainingEvents : integer := 5;
        variable frame : integer := 0;
        variable Li, Lo : line;
        variable itest, iobj : integer;
        variable part : particle;
    begin
        rst <= '1';
        wait for 5 ns;
        rst <= '0';
        start <= '0';
        wait until rising_edge(clk);
        while remainingEvents > 0 loop
            if not endfile(Fi) then
                readline(Fi, Li);
                read(Li, itest);
                read(Li, iobj); if (iobj > 0) then newevent <= '1'; else newevent <= '0'; end if;
                for i in 0 to NTKSECTORS*NTKFIBERS-1  loop
                    read(Li, iobj); part.pt   := (to_signed(iobj, 16));
                    read(Li, iobj); part.eta  := (to_signed(iobj, 10));
                    read(Li, iobj); part.phi  := (to_signed(iobj, 10));
                                    part.rest := (others => '0');
                    p_in(i) <= particle_to_w64(part);
                end loop;
                start <= '1';
             else
                remainingEvents := remainingEvents - 1;
                newevent <= '0';
                p_in <= (others => (others => '0'));
                start <= '1';
            end if;
           -- ready to dispatch ---
            wait until rising_edge(clk);
            -- write out the output --
            write(Lo, frame, field=>5);  
            write(Lo, string'(" 1 ")); 
            write(Lo, newevent_out); 
            write(Lo, string'(" ")); 
            for i in 0 to NREGIONS-1 loop
                part := w64_to_particle(p_out(i));
                write(Lo, to_integer(part.pt),   field => 5); 
                write(Lo, to_integer(part.eta),  field => 5); 
                write(Lo, to_integer(part.phi),  field => 5); 
            end loop;
            writeline(Fo, Lo);
            frame := frame + 1;
            --if frame >= 50 then finish(0); end if;
        end loop;
        wait for 50 ns;
        finish(0);
    end process;

    
end Behavioral;
