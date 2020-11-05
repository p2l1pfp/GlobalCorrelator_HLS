library std;
use std.textio.all;
use std.env.all;
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.std_logic_textio.all;

use work.regionizer_data.all;
use work.pattern_textio.all;


entity testbench is
--  Port ( );
end testbench;

architecture Behavioral of testbench is
    constant NPATTERNS_IN  : natural := NTKSECTORS*NTKFIBERS + NCALOSECTORS*NCALOFIBERS + NMUFIBERS;
    constant NPATTERNS_OUT : natural := NTKSORTED + NCALOSORTED + NMUSORTED;

    signal clk : std_logic := '0';
    signal rst : std_logic := '0';
    signal start, ready, idle, done : std_logic;
    signal newevent, newevent_out : std_logic;

    signal tk_in:  w64s(NTKSECTORS*NTKFIBERS-1 downto 0) := (others => (others => '0'));
    signal tk_out: w64s(NTKSORTED-1 downto 0) := (others => (others => '0'));
    signal calo_in:  w64s(NCALOSECTORS*NCALOFIBERS-1 downto 0) := (others => (others => '0'));
    signal calo_out: w64s(NCALOSORTED-1 downto 0) := (others => (others => '0'));
    signal mu_in:  w64s(NMUFIBERS-1 downto 0) := (others => (others => '0'));
    signal mu_out: w64s(NMUSORTED-1 downto 0) := (others => (others => '0'));

    file Fi : text open read_mode  is "input-emp.txt";
    file Fo : text open write_mode is "output-emp-vhdl_tb.txt";

begin
    clk  <= not clk after 1.25 ns;
    
    uut : entity work.full_regionizer_mux
        generic map(MU_ETA_CENTER => 460)
        port map(ap_clk => clk, 
                 ap_rst => rst, 
                 ap_start => start,
                 ap_ready => ready,
                 ap_idle =>  idle,
                 ap_done => done,
                 tracks_start => start,
                 tracks_newevent => newevent,
                 tracks_in_0_0_V => tk_in( 0),
                 tracks_in_0_1_V => tk_in( 1),
                 tracks_in_1_0_V => tk_in( 2),
                 tracks_in_1_1_V => tk_in( 3), 
                 tracks_in_2_0_V => tk_in( 4),
                 tracks_in_2_1_V => tk_in( 5),
                 tracks_in_3_0_V => tk_in( 6),
                 tracks_in_3_1_V => tk_in( 7),
                 tracks_in_4_0_V => tk_in( 8),
                 tracks_in_4_1_V => tk_in( 9), 
                 tracks_in_5_0_V => tk_in(10),
                 tracks_in_5_1_V => tk_in(11),
                 tracks_in_6_0_V => tk_in(12),
                 tracks_in_6_1_V => tk_in(13),
                 tracks_in_7_0_V => tk_in(14),
                 tracks_in_7_1_V => tk_in(15), 
                 tracks_in_8_0_V => tk_in(16),
                 tracks_in_8_1_V => tk_in(17),
                 calo_start => start,
                 calo_newevent => newevent,
                 calo_in_0_0_V => calo_in( 0),
                 calo_in_0_1_V => calo_in( 1),
                 calo_in_0_2_V => calo_in( 2),
                 calo_in_0_3_V => calo_in( 3), 
                 calo_in_1_0_V => calo_in( 4),
                 calo_in_1_1_V => calo_in( 5),
                 calo_in_1_2_V => calo_in( 6),
                 calo_in_1_3_V => calo_in( 7),
                 calo_in_2_0_V => calo_in( 8),
                 calo_in_2_1_V => calo_in( 9), 
                 calo_in_2_2_V => calo_in(10),
                 calo_in_2_3_V => calo_in(11),
                 mu_start => start,
                 mu_newevent => newevent,
                 mu_in_0_V => mu_in(0),
                 mu_in_1_V => mu_in(1),
                 tracks_out => tk_out,
                 calo_out   => calo_out,
                 mu_out     => mu_out,
                 newevent_out => newevent_out
             );
   

    runit : process 
        variable remainingEvents : integer := 5;
        variable patterns_in  : w64s(NPATTERNS_IN  - 1 downto 0);
        variable patterns_out : w64s(NPATTERNS_OUT - 1 downto 0);
        variable patterns_in_valid : std_logic;
        variable patterns_in_valid_old : std_logic := '0';
        variable patterns_out_valid : std_logic := '1';
        variable frame : integer := 0;
        variable Li, Lo : line;
        variable itest, iobj : integer;
        variable part : particle;
        variable gpart : glbparticle;
    begin
        rst <= '1';
        wait for 5 ns;
        rst <= '0';
        start <= '0';
        tk_in <= (others => (others => '0'));
        calo_in <= (others => (others => '0'));
        mu_in <= (others => (others => '0'));
        wait until rising_edge(clk);
        while remainingEvents > 0 loop
            start <= '1';
            if not endfile(Fi) then
                read_pattern_frame(FI, patterns_in, patterns_in_valid);
            else
                patterns_in_valid := '0';
                patterns_in := (others => (others => '0'));
                remainingEvents := remainingEvents - 1;
            end if;
            if patterns_in_valid = '1' and patterns_in_valid_old = '0' then
                newevent <= '1';
            else
                newevent <= '0';
            end if;
            patterns_in_valid_old := patterns_in_valid;
            for i in 0 to NTKSECTORS*NTKFIBERS-1 loop
                tk_in(i) <= patterns_in(i);
            end loop;
            for i in 0 to NCALOSECTORS*NCALOFIBERS-1 loop
                calo_in(i) <= patterns_in(i+NTKSECTORS*NTKFIBERS);
            end loop;
            for i in 0 to NMUFIBERS-1 loop
                mu_in(i) <= patterns_in(i+NTKSECTORS*NTKFIBERS+NCALOSECTORS*NCALOFIBERS);
            end loop;
           -- ready to dispatch ---
            wait until rising_edge(clk);
            -- write out the output --
            write(Lo, frame, field=>5);  
            for i in 0 to NTKSORTED-1 loop
                patterns_out(i) := tk_out(i);
            end loop;
            for i in 0 to NCALOSORTED-1 loop
                patterns_out(i+NTKSORTED) := calo_out(i);
            end loop;
            for i in 0 to NMUSORTED-1 loop
                patterns_out(i+NTKSORTED+NCALOSORTED) := mu_out(i);
            end loop;
            write_pattern_frame(Fo, frame, patterns_out, patterns_out_valid);
            --if frame >= 50 then finish(0); end if;
            frame := frame + 1;
        end loop;
        wait for 50 ns;
        finish(0);
    end process;

    
end Behavioral;
